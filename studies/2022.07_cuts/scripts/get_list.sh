# juliet stuff
# ls /home/mmeier/data/simulations/table_based_sim/juliet/nue/high_energy/l2_prime/1/*.zst > files_nue_high_energy.txt
# ls /home/mmeier/data/simulations/table_based_sim/juliet/numu/high_energy/l2_prime/1/*.zst > files_numu_high_energy.txt
# ls /home/mmeier/data/simulations/table_based_sim/juliet/nutau/high_energy/l2_prime/1/*.zst > files_nutau_high_energy.txt
# ls /home/mmeier/data/simulations/table_based_sim/juliet/mu/high_energy/l2_prime/1/*.zst > files_mu_high_energy.txt
# ls /home/mmeier/data/simulations/table_based_sim/juliet/tau/high_energy/l2_prime/1/*.zst > files_tau_high_energy.txt

# ls /home/mmeier/data/simulations/table_based_sim/juliet/nue/very_high_energy/l2_prime_cteq/1/*.zst > files_nue_high_energy.txt
# ls /home/mmeier/data/simulations/table_based_sim/juliet/numu/very_high_energy/l2_prime_cteq/1/*.zst > files_numu_high_energy.txt
# ls /home/mmeier/data/simulations/table_based_sim/juliet/nutau/very_high_energy/l2_prime_cteq/1/*.zst > files_nutau_high_energy.txt
# ls /home/mmeier/data/simulations/table_based_sim/juliet/mu/very_high_energy/l2_prime_cteq/1/*.zst > files_mu_high_energy.txt
# ls /home/mmeier/data/simulations/table_based_sim/juliet/tau/very_high_energy/l2_prime_cteq/1/*.zst > files_tau_high_energy.txt

# hdfwriter-merge /disk20/users/brian/IceCube/juliet/nue_high_energy/1/*.hdf5 -o /disk20/users/brian/IceCube/juliet/nue_high_energy_merged_998files.hdf5
# hdfwriter-merge /disk20/users/brian/IceCube/juliet/numu_high_energy/1/*.hdf5 -o /disk20/users/brian/IceCube/juliet/numu_high_energy_merged_999files.hdf5
# hdfwriter-merge /disk20/users/brian/IceCube/juliet/nutau_high_energy/1/*.hdf5 -o /disk20/users/brian/IceCube/juliet/nutau_high_energy_merged_996files.hdf5
# hdfwriter-merge /disk20/users/brian/IceCube/juliet/mu_high_energy/1/*.hdf5 -o /disk20/users/brian/IceCube/juliet/mu_high_energy_merged_998files.hdf5
# hdfwriter-merge /disk20/users/brian/IceCube/juliet/tau_high_energy/1/*.hdf5 -o /disk20/users/brian/IceCube/juliet/tau_high_energy_merged_999files.hdf5

# nugen stuff

# all files
# ls /disk20/users/brian/IceCube/nugen/21218/*/*.zst > files_21218.txt
# ls /disk20/users/brian/IceCube/nugen/21220/*/*.zst > files_21220.txt
# ls /disk20/users/brian/IceCube/nugen/21221/*/*.zst > files_21221.txt

# just a subrun
# ls /disk20/users/brian/IceCube/nugen/21218/0000000-0000999/*.zst > files_21218.txt
# ls /disk20/users/brian/IceCube/nugen/21220/0000000-0000999/*.zst > files_21220.txt
# ls /disk20/users/brian/IceCube/nugen/21221/0000000-0000999/*.zst > files_21221.txt

# hdfwriter-merge /disk20/users/brian/IceCube/nugen/hdf5/21218/*/*.hdf5 -o /disk20/users/brian/IceCube/nugen/hdf5/merged_21218_11991files.hdf5
# hdfwriter-merge /disk20/users/brian/IceCube/nugen/hdf5/21220/*/*.hdf5 -o /disk20/users/brian/IceCube/nugen/hdf5/merged_21220_10000files.hdf5
# hdfwriter-merge /disk20/users/brian/IceCube/nugen/hdf5/21221/*/*.hdf5 -o /disk20/users/brian/IceCube/nugen/hdf5/merged_21221_10000files.hdf5


# corsika stuff

# ls /disk20/users/brian/IceCube/corsika/level4/21962/0000000-0000999/* > files_21962.txt
# ls /disk20/users/brian/IceCube/corsika/level4/22023/0000000-0000999/* > files_22023.txt

#hdfwriter-merge /disk20/users/brian/IceCube/corsika/level4/hdf5/21962/*/*_p.hdf5 -o /disk20/users/brian/IceCube/corsika/level4/hdf5/corsika_21962_merged_p_998files.hdf5
#hdfwriter-merge /disk20/users/brian/IceCube/corsika/level4/hdf5/22023/*/*_p.hdf5 -o /disk20/users/brian/IceCube/corsika/level4/hdf5/corsika_22023_merged_p_999files.hdf5
#hdfwriter-merge /disk20/users/brian/IceCube/corsika/level4/hdf5/21962/*/*_fe.hdf5 -o /disk20/users/brian/IceCube/corsika/level4/hdf5/corsika_21962_merged_fe_998files.hdf5
#hdfwriter-merge /disk20/users/brian/IceCube/corsika/level4/hdf5/22023/*/*_fe.hdf5 -o /disk20/users/brian/IceCube/corsika/level4/hdf5/corsika_22023_merged_fe_999files.hdf5
