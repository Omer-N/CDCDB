BASE_DIR=""
now=$(date +'%d.%m.%Y')
mkdir -p data/final_schema/"$now"

aact_with_identifiers_path="data/final_schema/${now}/aact_combs__with_identifiers.csv"
echo 'Creating new version for C-DCDB'
echo 'Current Date' "$now"

if python "$BASE_DIR"clinical_trials_combinations.py --aact_params_file_path input_data/aact_credentials.json "data/final_schema/${now}/aact_combs.csv"; then
  echo 'Create clinical_trials_combinations'
else
  echo 'Failed creating clinical_trials_combinations'
  exit 1
fi
#
echo 'Adding identifiers for drugs for AACT'
if python add_identifier_to_df.py --as_str_array True "data/final_schema/${now}/aact_combs.csv" selected_name identifiers_entity "$aact_with_identifiers_path"; then
  echo 'Successfully added identifiers'
else
  echo 'Failed adding identifiers for AACT'
  exit 1
fi

echo 'Transforming data (normalizing)'
if python schema_transforming.py "$aact_with_identifiers_path" data/final_schema/"${now}" input_data/dbid_to_compound.csv; then
  echo 'Successfully transformed (normalized) the data'
else
  echo 'Failed to transform (normalize) the data'
  exit 1
fi

if python orange_book.py input_data/orange_book "data/final_schema/${now}/orangebook_combs_df.csv"; then
  echo 'Created orangebook tables successfully'
else
  echo 'Failed to create orangebook tables'
  exit 1
fi

cp input_data/patents/* data/final_schema/"${now}"

if python create_unnormalized_combs_db.py data/final_schema/"${now}" data/final_schema/"${now}"; then
  echo 'Created unnormalized version successfully'
else
  echo 'Failed to create unnormalized version'
  exit 1
fi

(
  cd data/final_schema/"${now}"
  zip -r "${now}.zip" .
)

cp "data/final_schema/${now}/${now}.zip" ../../../versions

cp data/final_schema/"${now}"/web_preview.csv ../../../latestVersion
echo $(date +'%d %b %Y') >../../../latestVersion/date.txt
