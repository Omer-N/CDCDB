import json
import pandas as pd
from tqdm import tqdm
import sys
sys.path.insert(0, '../..')

from src.drug_combs.add_identifier_to_df import DataframeDrugIdentifiersAdder
from src.drug_identfiers_resolver.identifiers_resolver import DrugIdentifiersResolver, APIBasedIdentifiersResolver, \
    WikiDataIdsResolver


class OrangeBookParser(object):
    def __init__(self, orange_book_path):
        self.orange_book_path = orange_book_path

    def get_combs(self) -> pd.DataFrame:
        results = []
        appl_no_to_patent_info = self.get_appl_no_to_product()

        with open(f'{str(self.orange_book_path)}/products.txt') as orange_book_file:
            for orange_book_line in orange_book_file.readlines():
                splitted = orange_book_line.split("~")
                ingredients = splitted[0]
                trade_name = splitted[2]
                appl_type = splitted[5]
                appl_no = splitted[6]
                product_no = splitted[7]
                te_code = splitted[8]
                approval_date = splitted[9]
                rld = splitted[10]
                rs = splitted[11]
                orangebook_type = splitted[12]
                applicant_full_name = splitted[13].strip()

                if ';' in ingredients:
                    drugs = ingredients.split("; ")
                    combs = {'all_drugs_in_combination': drugs,
                             'Trade_Name': trade_name,
                             'Appl_Type': appl_type,
                             'Product_no': product_no,
                             'TE_Code': te_code,
                             'Approval_Date': approval_date,
                             'RLD': rld,
                             'RS': rs,
                             'TYPE': orangebook_type,
                             'Applicant_Full_Name': applicant_full_name,
                             }
                    patent_data = appl_no_to_patent_info.get(appl_no, None)
                    if patent_data is not None:
                        for k, v in patent_data.items():
                            combs[k] = v
                    for idx, drug in enumerate(drugs):
                        combs[f'drug_{idx}'] = drug
                    results.append(combs)
        result_df = pd.DataFrame(results[1:])
        result_df['all_drugs_in_combination'] = result_df['all_drugs_in_combination'].apply(json.dumps)
        result_df = result_df.drop_duplicates(['all_drugs_in_combination'])
        result_df['all_drugs_in_combination'] = result_df['all_drugs_in_combination'].apply(json.loads)
        return result_df

    def get_appl_no_to_product(self):
        appl_no_to_patent_info = {}
        with open(f'{str(self.orange_book_path)}/patent.txt') as orange_book_file:
            for orange_book_line in orange_book_file.readlines():
                splitted = orange_book_line.split("~")
                appl_no = splitted[1]
                product_no = splitted[2]
                patent_no = splitted[3]
                patent_expire_date = splitted[4]
                drug_substance_flag = splitted[5]
                drug_product_flag = splitted[6]
                patent_use_code = splitted[7]
                delist_flag = splitted[8]
                patent_submitted_date = splitted[9]
                appl_no_to_patent_info[appl_no] = {"product_no": product_no, 'patent_no': patent_no,
                                                   'patent_expire_date': patent_expire_date,
                                                   'substance_flag': drug_substance_flag,
                                                   'product_flag': drug_product_flag,
                                                   'patent_use_code': patent_use_code,
                                                   'delist_flag': delist_flag,
                                                   'patent_submitted_date': patent_submitted_date}
        return appl_no_to_patent_info


def get_drugs_identifiers_adder():
    wikidata_ids_resolver = WikiDataIdsResolver("../drug_combs/input_data/qid_to_drugbank.json",
                                                "../drug_combs/input_data/qid_to_pubchem.json",
                                                cache_file_path="wikidata_disk_cache")
    resolver = DrugIdentifiersResolver(wikidata_ids_resolver, APIBasedIdentifiersResolver())
    drug_identifiers_adder = DataframeDrugIdentifiersAdder(resolver, False)
    return drug_identifiers_adder


def main(args):
    orange_book_parser = OrangeBookParser(args.orangebook_path)
    combs_df = orange_book_parser.get_combs()
    combs_df = combs_df.fillna("")
    drug_identifiers_adder = get_drugs_identifiers_adder()
    for column in combs_df.columns:
        if column.startswith("drug_"):
            combs_df = drug_identifiers_adder.add_identifiers_column(combs_df, column, f'identifiers_{column}')

    columns = combs_df.columns
    result_df = []
    for row in tqdm(combs_df.iterrows(), desc="Removing empty identifiers: rows"):
        row = row[1]
        drugs_names = []
        pbids = []
        dbids = []
        for column in columns:
            if column.startswith("drug_"):
                if row[column] == '':
                    row[f"identifiers_{column}"] = ''
                    continue
                drug_i_name = str(row[column]).strip()
                drugs_names.append(drug_i_name)
            elif column.startswith("identifiers_"):
                if row[column] == '':
                    continue
                dbids.append(row[column][0])
                pbids.append(row[column][1])
        row['drugs_names'] = json.dumps(drugs_names)
        row['drugbank_ids'] = json.dumps(dbids)
        row['pubchem_ids'] = json.dumps(pbids)
        result_df.append(row)

    result_df = pd.DataFrame(result_df)
    cols_to_drop = [x for x in result_df.columns if x.startswith("drug_")] + [x for x in result_df.columns if
                                                                              x.startswith("identifiers_")] + ['all_drugs_in_combination']
    combs_df = result_df.drop(cols_to_drop, axis=1)

    if args.parsed_df_path[-3:] == 'csv':
        combs_df.to_csv(f'{args.parsed_df_path}', index=False)
    elif args.parsed_df_path[-3:] == 'xlsx':
        combs_df.to_excel(f'{args.parsed_df_path}', index=False)
    else:
        exit("Unsupported export format, only csv, xlsx allowed")


if __name__ == '__main__':
    import argparse

    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("orangebook_path", type=str)
    argument_parser.add_argument("parsed_df_path", type=str, help="path to output the combs df [csv/xlsx]")
    args = argument_parser.parse_args()
    main(args)

