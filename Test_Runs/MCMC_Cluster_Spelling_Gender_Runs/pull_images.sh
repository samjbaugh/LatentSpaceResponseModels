seeds=$(cat seed_list.txt)
x=0
for seednum in $seeds
do
	x=$(($x+1))
	#cp Images/plots_config_1_seed_"$seednum"_data_spellingGender/initial_configuration.png Images_Consolidated/initial_config_run_"$x".png
	cp Images/plots_config_1_seed_"$seednum"_data_spellingGender/iteration_5000.png Images_Consolidated/iteration_5000_run_"$x".png
done




