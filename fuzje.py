import re
from os.path import isfile
import csv
import copy
from statistics import mean

import wget
import ast


# be careful - many function in this script require csv file with specific column names and file names

# creates csv with knot cores basing on data from https://alphaknot.cent.uw.edu.pl/
def csv_with_knot_cores(family_id):
    with open(f"Kopia fuzje-results - {family_id}.csv") as file_csv:
        rows = list(csv.DictReader(file_csv))
        csv_columns = list(rows[0].keys())
        csv_columns.append("Knot core - Alphafold")
        csv_columns.append("Knot core - Rosetta")
        for protein in rows:
            id_protein = protein["Uniprot ID"]
            protein["Knot core - Alphafold"] = ""
            protein["Knot core - Rosetta"] = ""
            if not isfile(f"{id_protein}_a.html"):
                downloading_html(protein["AlphaFold"], "a", id_protein)
            if not isfile(f"{id_protein}_r.html"):
                downloading_html(protein["Rosetta"], "r", id_protein)
            if isfile(f"{id_protein}_a.html"):
                protein["Knot core - Alphafold"] = knot_core_from_html_file_alphaknot(f"{id_protein}_a.html")
            if isfile(f"{id_protein}_r.html"):
                protein["Knot core - Rosetta"] = knot_core_from_html_file_alphaknot(f"{id_protein}_r.html")
        with open(f"{family_id}.csv", 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
            writer.writeheader()
            for data in rows:
                writer.writerow(data)


# finds data about knot cores in html file (site https://alphaknot.cent.uw.edu.pl/)
def knot_core_from_html_file_alphaknot(file_name):
    with open(file_name) as file:
        markers = re.findall(r'(?<=<div id=marker_)[\w#-]+[^></div]', file.read())
        core = ""
        for marker in markers:
            core += marker + "\n"
        return core[:-1]


# download html file
def downloading_html(string_with_link_in_it, r_or_a, name):
    try:
        wget.download(re.search(r'(https?)[\SA-Za-z/d]*[^"\s]', string_with_link_in_it).group(0),
                      f"{name}_{r_or_a}.html")
    except:
        pass


# downloads matrix of knots from https://alphaknot.cent.uw.edu.pl/
def downloading_matrixes(link_alphaknot, matrix_name):
    compute = link_alphaknot.find("compute")
    link = link_alphaknot[:compute] + "compute_static" + link_alphaknot[compute + len("compute"):-1] * 2
    link2 = link_alphaknot[:compute] + "compute_static" + link_alphaknot[compute + len("compute"):-1]
    if isfile(matrix_name):
        return True
    try:
        wget.download(f"{link}_1.txt", matrix_name)
        return True
    except:
        try:
            wget.download(f"{link2}_1.txt", matrix_name)
            return True
        except:
            return False


# reading knot matrix file to find knot cores
def getting_knot_core_from_matrix(matrix_file_name):
    with open(matrix_file_name) as matrix_file:
        data = matrix_file.read()
        matrix = ast.literal_eval(data)
        result = {}
        for subchain, knots in matrix.items():
            for knot_name, probability in knots.items():
                if probability >= 0.5 and knot_name != "0_1" and knot_name != "3_1#3_1|8_20":
                    if knot_name not in result.keys():
                        pom = set()
                        pom.add((subchain, probability))
                    else:
                        pom = result[knot_name].copy()
                        for el in result[knot_name]:
                            start = subchain[0]
                            end = subchain[1]
                            '''if start > el[0][1]:
                                pom.add((subchain, probability))
                            elif end < el[0][0]:
                                pom.add((subchain, probability))'''
                            if start >= el[0][0] and end <= el[0][1]:
                                pom.remove(el)
                            pom.add((subchain, probability))
                            '''else:
                                pom.add((subchain, probability))'''
                    result[knot_name] = pom
    answer = copy.deepcopy(result)
    for knot_type, set_of_subchains in result.items():
        for el in set_of_subchains:
            for el2 in set_of_subchains:
                if el[0][0] in range(el2[0][0], el2[0][1] + 1) or el[0][1] in range(el2[0][0], el2[0][1] + 1):
                    if el[0][1] - el[0][0] < el2[0][1] - el2[0][0]:
                        answer[knot_type].discard(el2)
                    if el[0][1] - el[0][0] == el2[0][1] - el2[0][0]:
                        if el[1] > el2[1]:
                            answer[knot_type].discard(el2)
                        if el[1] < el2[1]:
                            answer[knot_type].discard(el)
    text = ""
    for knot_type, subchains in answer.items():
        text += f'{knot_type}: '
        for subchain, probability in subchains:
            text += f'{subchain[0]} - {subchain[1]}: {probability}; '
        text += "\n"
    return text[:-1]


# reading knot matrix file to find knot cores and return mean of start and end positions
def getting_mean_of_knot_cores_from_matrix(matrix_file_name):
    with open(matrix_file_name) as matrix_file:
        data = matrix_file.read()
        matrix = ast.literal_eval(data)
        result = {}
        for subchain, knots in matrix.items():
            for knot_name, probability in knots.items():
                # there will be cases with very high probability but also with to many unimportant aminoacids
                if 0.5 <= probability <= 0.6 and knot_name != "0_1" and knot_name != "3_1#3_1|8_20":
                    if knot_name not in result.keys():
                        result[knot_name] = [[subchain]]
                    else:
                        pom = copy.deepcopy(result)
                        for i in range(0, len(result[knot_name])):
                            start = subchain[0]
                            end = subchain[1]
                            starts = [j[0] for j in pom[knot_name][i]]
                            ends = [j[1] for j in pom[knot_name][i]]
                            if start > max(ends) or end < min(starts):
                                # this case means that we are dealing with another knot of this type
                                pom[knot_name].append([subchain])
                            else:
                                pom[knot_name][i].append(subchain)
                        result = pom
    text = ""
    for knot_type, ranges in result.items():
        text += f"{knot_type}: "
        for el in ranges:
            starts = [i[0] for i in el]
            ends = [i[1] for i in el]
            text += f'{round(mean(starts))} - {round(mean(ends))} '
            text += "\n"
    return text[:-1]


# makes csv file with new columns with knot cores based on matrix from https://alphaknot.cent.uw.edu.pl/
def knot_cores_in_family_creates_csv(family_id):
    with open(f"{family_id}.csv") as file_csv:
        rows = list(csv.DictReader(file_csv))
        col_name = "Alphafold knot cores if many knot types"
        for protein in rows:
            id = protein["Uniprot ID"]
            knots = protein["Knot core - Alphafold"].split()
            if len(knots) > 1 or len(protein["Knot core - Alphafold"].split("#")) > 1:
                link = re.search(r'(https?)[\SA-Za-z/d]*[^"\s]', protein["AlphaFold"]).group(0)
                if downloading_matrixes(link, f"{id}_matrix.txt"):
                    protein[col_name] = getting_knot_core_from_matrix(f"{id}_matrix.txt")
                else:
                    protein[col_name] = "problem with getting matrix"
            else:
                protein[col_name] = ""
        with open(f"{family_id}_with_knot_cores.csv", 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            for row in rows:
                writer.writerow(row)


# makes csv file with new columns with knot cores means based on matrix from https://alphaknot.cent.uw.edu.pl/
def knot_cores_mean_in_family_creates_csv(family_id):
    with open(f"{family_id}_domains.csv") as file_csv:
        rows = list(csv.DictReader(file_csv))
        col_name = "Alphafold knot cores - mean"
        for protein in rows:
            id = protein["Uniprot ID"]
            knots = protein["Knot core - Alphafold"].split()
            if len(knots) > 1 or len(protein["Knot core - Alphafold"].split("#")) > 1:
                link = re.search(r'(https?)[\SA-Za-z/d]*[^"\s]', protein["AlphaFold"]).group(0)
                if downloading_matrixes(link, f"{id}_matrix.txt"):
                    protein[col_name] = getting_mean_of_knot_cores_from_matrix(f"{id}_matrix.txt")
                else:
                    protein[col_name] = "problem with getting matrix"
            else:
                protein[col_name] = ""
        with open(f"{family_id}_with_knot_cores_means.csv", 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            for row in rows:
                writer.writerow(row)


# finding domain ranges in proteins from the family (informations about domains based on PFAM database)
def getting_domain_range_from_pfam(protein_id):
    if not isfile(f"{protein_id}_pfam.html"):
        try:
            wget.download(f"http://pfam.xfam.org/protein/{protein_id}#tabview=tab0", f"{protein_id}_pfam.html")
        except:
            return "no domains in pfam\n"
    domains = ""
    with open(f"{protein_id}_pfam.html") as file:
        lines = file.readlines()
        for i in range(0, len(lines)):
            domain = re.search(r'(?<=<td class="pfama_)[A-Za-z\d]*', lines[i])
            if domain:
                i = i + 2
                start = re.search(r'\d+', lines[i]).group(0)
                i = i + 1
                end = re.search(r'\d+', lines[i]).group(0)
                domains += f"{domain.group(0)}: {start} - {end}\n"
    return domains


# creates new csv file with added column with domain ranges
def domain_to_id_in_family(family_id):
    with open(f"{family_id}_with_knot_cores.csv") as file_csv:
        rows = list(csv.DictReader(file_csv))
        col_name = "domains ranges"
        for protein in rows:
            id = protein["Uniprot ID"]
            protein[col_name] = getting_domain_range_from_pfam(id)[:-1]
        with open(f"{family_id}_domains.csv", 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            for row in rows:
                writer.writerow(row)


def main():
    list_of_families = ["PF01699-PF01699_PF01699-PF01699", "PF00588_PF00588", "PF03587_PF03587"]
    for family in list_of_families:
        '''csv_with_knot_cores(family)
        knot_cores_in_family_creates_csv(family)
        domain_to_id_in_family(family)'''
        knot_cores_mean_in_family_creates_csv(family)


main()
# knot_cores_in_family_creates_csv("PF00588_PF00588")
