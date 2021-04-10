import pandas as pd
from os import sep


def main(args):
    orangebook_df = pd.read_csv(f'{args.input_dir}{sep}orangebook_combs_df.csv',
                                usecols=['drugs_names', 'drugbank_ids', 'pubchem_ids'])
    orangebook_df['source_id'] = 'N/A'
    orangebook_df['source'] = 'orangebook'
    orangebook_df.columns = ['drugs', 'drugbank_identifiers', 'pubchem_identifiers', 'source_id', 'source']
    aact_df = pd.read_csv(f'{args.input_dir}{sep}design_group_df.csv', usecols=['selected_name',
                                                                                'drugbank_identifier',
                                                                                'pubchem_identifier', 'nct_id'])
    aact_df['source_id'] = aact_df['nct_id']
    aact_df = aact_df.drop('nct_id', axis=1)
    aact_df['source'] = 'clinicaltrials.gov'
    aact_df.columns = ['drugs', 'drugbank_identifiers', 'pubchem_identifiers', 'source_id', 'source']

    patents_df = pd.read_csv(f'{args.input_dir}{sep}transformed_patents_drug.csv',
                             usecols=['drugs_names', 'drugbank_identifiers', 'pubchem_identifiers', 'Patent ID'])
    patents_df['source_id'] = patents_df['Patent ID']
    patents_df = patents_df.drop('Patent ID', axis=1)
    patents_df['source'] = 'patents'
    patents_df.columns = ['drugs', 'drugbank_identifiers', 'pubchem_identifiers', 'source_id', 'source']

    all_combs = pd.concat([aact_df, patents_df, orangebook_df])
    all_combs.to_csv(f'{args.output_path}/all_combs_unormalized.csv', index=False)

    create_web_preview_table(all_combs, f'{args.output_path}/web_preview.csv')


def create_web_preview_table(all_combs, web_preview_path):
    web_preview = all_combs
    web_preview['drugbank_identifiers'] = web_preview['drugbank_identifiers'].apply(lambda x: eval(x))
    web_preview['pubchem_identifiers'] = web_preview['pubchem_identifiers'].apply(lambda x: eval(x))
    web_preview['drugs'] = web_preview['drugs'].apply(lambda x: eval(x))

    def add_best_match_name(row):
        result = []
        dbid_identifiers = row['drugbank_identifiers']
        names = row['drugs']
        for i in range(len(dbid_identifiers)):
            current_dbid = dbid_identifiers[i]
            if current_dbid != '-1' and current_dbid in dbid_to_name:
                result.append(dbid_to_name[current_dbid])
            else:
                result.append(names[i][0])
        return ','.join(result)

    # TODO: add drugbank name
    drugbank_names_df = pd.read_csv('input_data/drugbank_drug_names.csv')
    dbid_to_name = drugbank_names_df.set_index('drugBank_id').to_dict()['Drug name']
    web_preview['drugs'] = web_preview.apply(lambda row: add_best_match_name(row), axis=1)

    def transform_identifiers_for_web(arr):
        arr = list(map(lambda x: x if x != '-1' else 'NA', arr))
        return ';'.join(arr)

    web_preview['drugbank_identifiers'] = web_preview['drugbank_identifiers'].apply(
        lambda x: transform_identifiers_for_web(x))
    web_preview['pubchem_identifiers'] = web_preview['pubchem_identifiers'].apply(
        lambda x: transform_identifiers_for_web(x))

    web_preview.drop_duplicates().to_csv(web_preview_path, index=False)


if __name__ == '__main__':
    import argparse

    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("input_dir", help="input dir")
    argument_parser.add_argument("output_path", help="output_file.csv")
    args = argument_parser.parse_args()
    main(args)
