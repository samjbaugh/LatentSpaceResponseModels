mv MCMC_Template/Images/* ../OutputArchiveLatentSpace/Images
mv MCMC_Template/Saved_output/* ../OutputArchiveLatentSpace/Saved_output
for x in $(ls Test_Runs)
do
	mv Test_Runs/$x/Images/* ../OutputArchiveLatentSpace/Images
	mv Test_Runs/$x/Saved_output/* ../OutputArchiveLatentSpace/Saved_output
done


