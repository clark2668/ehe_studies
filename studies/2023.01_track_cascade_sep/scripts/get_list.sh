# juliet stuff


# juliet HE
# ls /home/mmeier/data/simulations/table_based_sim/juliet/numu/high_energy/L2/1/*.zst > files_numu_high_energy.txt
# ls /home/mmeier/data/simulations/table_based_sim/juliet/nutau/high_energy/L2/1/*.zst > files_nutau_high_energy.txt
# ls /home/mmeier/data/simulations/table_based_sim/juliet/tau/high_energy/L2/1/*.zst > files_tau_high_energy.txt


# juliet VHE
# ls /home/mmeier/data/simulations/table_based_sim/juliet/numu/very_high_energy/l2_prime/1/*.zst > files_numu_very_high_energy.txt
# ls /home/mmeier/data/simulations/table_based_sim/juliet/nutau/very_high_energy/l2_prime/1/*.zst > files_nutau_very_high_energy.txt
# ls /home/mmeier/data/simulations/table_based_sim/juliet/tau/very_high_energy/l2_prime/1/*.zst > files_tau_very_high_energy.txt

## old L2? (Maybe max will tell me where they were moved to...)
# ls /disk19/users/mmeier/simulations/table_based_sim/juliet/nue/high_energy_bugged/l2_prime/1/*.zst > files_nue_high_energy.txt
# ls /home/mmeier/data/simulations/table_based_sim/juliet/mu/high_energy/L2/1/*.zst > files_mu_high_energy.txt
# ls /home/mmeier/data/simulations/table_based_sim/juliet/nue/very_high_energy/l2_prime/1/*.zst > files_nue_very_high_energy.txt
# ls /home/mmeier/data/simulations/table_based_sim/juliet/mu/very_high_energy/l2_prime/1/*.zst > files_mu_very_high_energy.txt

hdfwriter-merge /disk20/users/brian/IceCube/juliet/nue_high_energy/1/*.hdf5 -o /disk20/users/brian/IceCube/juliet/nue_high_energy/nue_high_energy_merged_998files.hdf5
hdfwriter-merge /disk20/users/brian/IceCube/juliet/nue_very_high_energy/1/*.hdf5 -o /disk20/users/brian/IceCube/juliet/nue_very_high_energy/nue_very_high_energy_merged_999files.hdf5
hdfwriter-merge /disk20/users/brian/IceCube/juliet/mu_high_energy/1/*.hdf5 -o /disk20/users/brian/IceCube/juliet/mu_high_energy/mu_high_energy_merged_999files.hdf5
hdfwriter-merge /disk20/users/brian/IceCube/juliet/mu_very_high_energy/1/*.hdf5 -o /disk20/users/brian/IceCube/juliet/mu_very_high_energy/mu_very_high_energy_merged_999files.hdf5


## juliet L4
#ls /home/mmeier/data/EHE/level4_v2/nue_high_energy/1/*.zst > files_nue_high_energy.txt
#ls /home/mmeier/data/EHE/level4_v2/nue_very_high_energy/1/*.zst > files_nue_very_high_energy.txt
#ls /home/mmeier/data/EHE/level4_v2/mu_high_energy_01/1/*.zst > files_mu_high_energy.txt
#ls /home/mmeier/data/EHE/level4_v2/mu_very_high_energy/1/*.zst > files_mu_very_high_energy.txt

# hdfwriter-merge /disk20/users/brian/IceCube/juliet/nue_high_energy/1/*.hdf5 -o /disk20/users/brian/IceCube/juliet/nue_high_energy/nue_high_energy_merged_998files.hdf5
# hdfwriter-merge /disk20/users/brian/IceCube/juliet/nue_very_high_energy/1/*.hdf5 -o /disk20/users/brian/IceCube/juliet/nue_very_high_energy/nue_very_high_energy_merged_999files.hdf5
# hdfwriter-merge /disk20/users/brian/IceCube/juliet/mu_high_energy/1/*.hdf5 -o /disk20/users/brian/IceCube/juliet/mu_high_energy/mu_high_energy_merged_999files.hdf5
# hdfwriter-merge /disk20/users/brian/IceCube/juliet/mu_very_high_energy/1/*.hdf5 -o /disk20/users/brian/IceCube/juliet/mu_very_high_energy/mu_very_high_energy_merged_999files.hdf5

# just merge Max's files (duh....)

# hdfwriter-merge /home/mmeier/data/EHE/level4_v2/nue_high_energy/1/*0000000*.hd5 -o /disk20/users/brian/IceCube/juliet/nue_high_energy/nue_high_energy_merged_998files.hdf5
# hdfwriter-merge /home/mmeier/data/EHE/level4_v2/nue_very_high_energy/1/*0000000*.hd5 -o /disk20/users/brian/IceCube/juliet/nue_very_high_energy/nue_very_high_energy_merged_999files.hdf5
# hdfwriter-merge /home/mmeier/data/EHE/level4_v2/mu_high_energy_01/1/*0000000*.hd5 -o /disk20/users/brian/IceCube/juliet/mu_high_energy/mu_high_energy_merged_999files.hdf5
# hdfwriter-merge /home/mmeier/data/EHE/level4_v2/mu_very_high_energy/1/*0000000*.hd5 -o /disk20/users/brian/IceCube/juliet/mu_very_high_energy/mu_very_high_energy_merged_999files.hdf5
