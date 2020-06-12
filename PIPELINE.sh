clear

echo "SARS-CoV-2 Indel analysis pipeline"
echo "Alexander G. Lucaci"
echo "last update: 06 08 2020"
echo ""

#
#echo "(Controller) Installing HyPhy"
#bash install_hyphy.sh

#
echo "(Controller) launching step 1"
#saves a 'wg_variation_daybyday_coloredbycountry.csv
#bash step_1_indel_analysis_pipeline_mod.sh

#
echo "(Controller) launching step 2"
#bash step_2_genebygene_mod.sh

#
echo "(Controller) launching for_observable.sh"
bash for_observable.sh

#
#for_observable_withgenomicpos.sh

