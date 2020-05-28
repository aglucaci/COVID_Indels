clear

echo "SARS-CoV-2 Indel analysis pipeline"
echo "Alexander G. Lucaci"
echo "last update: 05 27 2020"
echo ""

#
bash step_1_indel_analysis_pipeline_mod.sh

#
bash step_2_genebygene_mod.sh

#
bash for_observable.sh

