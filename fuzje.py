import os
import re
import csv
import copy
from statistics import mean
import wget
import ast
from os.path import isfile


# be careful - many function in this script require a csv file with specific columns and file names
# required columns:
# "Uniprot ID" with ids in uniprot standard
# "AlphaFold" and "Rosetta" which must contains link to on alphaknot database, but can also contain more data
# "Knot core - Alphafold" to run knot_cores_in_family_creates_csv(family_id) or knot_cores_mean_in_family_creates_csv(family_id).
# You can create a csv file with the column named "Knot core - Alphafold" with function  csv_with_knot_cores(family_id).


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
        # matrix file contains a dict with subchains as key and dict of knot types and probability as value
        matrix = ast.literal_eval(data)
        result = {}
        for subchain, knots in matrix.items():
            for knot_name, probability in knots.items():
                if probability >= 0.5 and knot_name != "0_1" and knot_name != "3_1#3_1|8_20":
                    if knot_name not in result.keys():
                        pom = set()
                        pom.add((subchain, probability))
                    else:
                        pom = copy.deepcopy(result)[knot_name]
                        for el in result[knot_name]:
                            start = subchain[0]
                            end = subchain[1]
                            # this case means that we can have a shorter knot core in this range, so we can replace it
                            if start >= el[0][0] and end <= el[0][1]:
                                pom.remove(el)
                            pom.add((subchain, probability))
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
                            starts = [j[0] for j in result[knot_name][i]]
                            ends = [j[1] for j in result[knot_name][i]]
                            if start > max(ends) or end < min(starts):
                                # this case means that we are dealing with another knot of this type
                                print("i co ")
                                pom[knot_name].append([subchain])
                            else:
                                pom[knot_name][i].append(subchain)
                        result = pom
    text = ""
    print(len(result["3_1"]))
    for knot_type, ranges in result.items():
        text += f"{knot_type}: "
        for el in ranges:
            starts = [i[0] for i in el]
            ends = [i[1] for i in el]
            text += f'{round(mean(starts))} - {round(mean(ends))} '
            text += "\n"
    return text[:-1]


def mean_of_knot_cores(matrix_file_name, reference_knot_cores):
    with open(matrix_file_name) as matrix_file:
        data = matrix_file.read()
        matrix = ast.literal_eval(data)
        result = {key: [[e[0]] for e in s] for key, s in reference_knot_cores.items()}
        for subchain, knots in matrix.items():
            for knot_name, probability in knots.items():
                if probability > 0.5 and knot_name != "0_1" and knot_name != "3_1#3_1|8_20":
                    start = subchain[0]
                    end = subchain[1]
                    for i in range(0, len(reference_knot_cores[knot_name])):
                        chain, probability_reference = reference_knot_cores[knot_name][i]
                        ref_start = range(int(chain[0]) - 20, int(chain[0]) + 20)
                        ref_end = range(int(chain[1]) - 20, int(chain[1]) + 20)
                        if start in ref_start and end in ref_end:
                            if float(probability_reference) - 0.02 <= probability <= float(
                                    probability_reference) + 0.02:
                                result[knot_name][i].append(subchain)

    text = ""
    for knot_type, ranges in result.items():
        text += f"{knot_type}: "
        for el in ranges:
            starts = [int(i[0]) for i in el]
            ends = [int(i[1]) for i in el]
            text += f'{round(mean(starts))} - {round(mean(ends))}; '
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
                    reference_knot_core = protein["Alphafold knot cores if many knot types"].split("; \n")
                    dict_of_knots = {}
                    for knot in reference_knot_core:
                        knot = knot.strip()
                        knot_name = knot[:knot.find(":")]
                        if knot_name != "0_1" and knot_name != "3_1#3_1|8_20":
                            knots = knot[knot.find(":") + 2:].split("; ")
                            dict_of_knots[knot_name] = []
                            for el in knots:
                                numbers = re.findall(r"[\d\.]+", el)
                                dict_of_knots[knot_name].append([(numbers[0], numbers[1]), numbers[2]])
                    protein[col_name] = mean_of_knot_cores(f"{id}_matrix.txt", dict_of_knots)
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
            return "no domains in PFAM. Probably not correct protein id\n", {}
    domains = ""
    domains_dict = {}
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
                if domain.group(0) not in domains_dict.keys():
                    domains_dict[domain.group(0)] = [(start, end)]
                else:
                    domains_dict[domain.group(0)].append((start, end))
    return domains, domains_dict


# creates new csv file with added column with domain ranges
def domain_to_id_in_family(family_id):
    with open(f"{family_id}_with_knot_cores_means.csv") as file_csv:
        rows = list(csv.DictReader(file_csv))
        col_name = "domains ranges"
        col_name_2 = "domains and knot core cover"
        for protein in rows:
            id = protein["Uniprot ID"]
            domain_range, dict_of_domains = getting_domain_range_from_pfam(id)
            protein[col_name] = domain_range[:-1]
            protein[col_name_2] = domain_and_knot_core_correlation(dict_of_domains,
                                                                   protein["Alphafold knot cores - mean"])
        with open(f"{family_id}_domains_and_domain_cover.csv", 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            for row in rows:
                writer.writerow(row)


def domain_and_knot_core_correlation(domain_ranges, knot_cores):
    if len(knot_cores) == 0 or domain_ranges == {}:
        return ""
    knot_cores = knot_cores + " \n"
    knot_cores = knot_cores.split("; \n")
    result = ""
    for knot_type_line in knot_cores:
        knot_name = knot_type_line[:knot_type_line.find(":")]
        knot_ranges = knot_type_line[len(knot_name) + 2:]
        knot_ranges = knot_ranges.split(";")
        for range_knot in knot_ranges[:-1]:
            start, end = re.findall(r"[\d\.]+", range_knot)
            knot_range = range(int(start), int(end) + 1)
            result += f"{knot_name}: "
            domains = list(domain_ranges.values())[0]
            for domain in domains:
                domain_range = range(int(domain[0]), int(domain[1]) + 1)
                covered = set(knot_range) & set(domain_range)
                if len(covered) > 0:
                    start_covered = min(covered)
                    end_covered = max(covered)
                    result += f"{start_covered} - {end_covered} (domena {domain[0]}, {domain[1]}), "
            result = result[:-2] + ";\n"
    print(result)
    return result[:-1]


def downloading_pdb(link_alphaknot, pdb_name, family_id):
    compute = link_alphaknot.find("compute")
    link = link_alphaknot[:compute] + "compute_static" + link_alphaknot[compute + len("compute"):-1] * 2
    link2 = link_alphaknot[:compute] + "compute_static" + link_alphaknot[compute + len("compute"):-1]
    directory = f"pdb_files_{family_id}"
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)
    if not isfile(pdb_name):
        try:
            wget.download(f"{link}_1.pdb", pdb_name)
        except:
            try:
                wget.download(f"{link2}_1.pdb", pdb_name)
            except:
                pass


def downloading_pdb_fol_all_proteins_in_family(family_id):
    with open(f"{family_id}.csv") as file_csv:
        rows = list(csv.DictReader(file_csv))
        for protein in rows:
            try:
                link = re.search(r'(https?)[\SA-Za-z/d]*[^"\s]', protein["AlphaFold"]).group(0)
                id = protein["Uniprot ID"]
                downloading_pdb(link, f"{id}.pdb", family_id)
            except:
                pass

def main():
    list_of_families = ["PF01699-PF01699_PF01699-PF01699", "PF00588_PF00588", "PF03587_PF03587"]
    for family in list_of_families:
        # csv_with_knot_cores(family)
        # knot_cores_in_family_creates_csv(family)
        # domain_to_id_in_family(family)
        downloading_pdb_fol_all_proteins_in_family(family)
        # knot_cores_mean_in_family_creates_csv(family)


main()
