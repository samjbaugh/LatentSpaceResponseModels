mv -r MCMC\ Template/Images/* ../OutputArchiveLatentSpace/Images
mv -r MCMC\ Template/Saved_Output/* ../OutputArchiveLatentSpace/SavedOutput
for x in $(ls Test\ Runs)
do
	mv -r Test\ Runs/$x/Images/* ../OutputArchiveLatentSpace/Images
	mv -r Test\ Runs/$x/Saved_Output/* ../OutputArchiveLatentSpace/SavedOutput
done


