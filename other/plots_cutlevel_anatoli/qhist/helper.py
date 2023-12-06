import os
from pypet import Trajectory
from neutrino_level.steps.general_modules.utils import add_dict_to_traj
from neutrino_level.steps.general_modules.label_maker import add_label_from_components_to_df
import numpy as np


def build_trajectory():
    traj = Trajectory(name='Trajectory', add_time=False)
    if os.uname().nodename in ['Brians-MBP', 'dhcp-172-16-174-147.private.pa.umd.edu']:
        DATA_PATH='/Users/brianclark/Documents/work/IceCube/ehe/new_i3files'
        JULIET_PATH = '/Users/brianclark/Documents/work/IceCube/ehe/new_i3files'
    elif os.uname().nodename in ['condor00']:
        DATA_PATH='/data/i3home/baclark/transfer_files'
        JULIET_PATH = '/data/i3home/baclark/transfer_files'
    else:
        raise NotImplementedError

    #####
    # numu high energy
    #####
    
    comp_dict_jnm = {}
    comp_dict_jnm['type'] = ['juliet_numu', 'Component type.']

    comp_dict_jnm['file_list'] = [
        [os.path.join(JULIET_PATH, 'numu_high_energy_merged_1kfiles.hdf5')],
        'Files for the component.']
    comp_dict_jnm['flux_names'] = [
        ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
         'ahlers_2010_1E18'],
        'Weight cols.']
    comp_dict_jnm['additional_obs'] = [['I3MCWeightDict.InteractionType',
                                       'NuMuProperties.muon_energy',
                                       'NuMuProperties.muon_energy_losses'],
                                      'Additional observables.']
    comp_dict_jnm['n_files_total'] = [10000, 'Number of files in Madison.']
    comp_dict_jnm['n_files_loaded'] = [[1000], 'Number of files loaded.']

    add_dict_to_traj(traj, 'components.numu_he', 'NuMu Component.',
                     comp_dict_jnm, to_log=False)


    #####
    # numu high energy
    #####

    comp_dict_jnmv = {}
    comp_dict_jnmv['type'] = ['juliet_numu', 'Component type.']

    comp_dict_jnmv['file_list'] = [
        [os.path.join(JULIET_PATH, 'numu_very_high_energy_merged_1kfiles.hdf5')],
        'Files for the component.']
    comp_dict_jnmv['flux_names'] = [
        ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
         'ahlers_2010_1E18'],
        'Weight cols.']
    comp_dict_jnmv['additional_obs'] = [['I3MCWeightDict.InteractionType',
                                       'NuMuProperties.muon_energy',
                                       'NuMuProperties.muon_energy_losses'],
                                      'Additional observables.']
    comp_dict_jnmv['n_files_total'] = [2000, 'Number of files in Madison.']
    comp_dict_jnmv['n_files_loaded'] = [[1000], 'Number of files loaded.']

    add_dict_to_traj(traj, 'components.numu_vhe', 'VHE NuMu Component.',
                     comp_dict_jnmv, to_log=False)
    
    #####
    # mu high energy
    #####
    
    comp_dict_jm = {}
    comp_dict_jm['type'] = ['juliet_mu', 'Component type.']

    comp_dict_jm['file_list'] = [
        [os.path.join(JULIET_PATH, 'mu_high_energy_merged_1kfiles.hdf5')],
        'Files for the component.']
    comp_dict_jm['flux_names'] = [
        ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
         'ahlers_2010_1E18'],
        'Weight cols.']
    comp_dict_jm['additional_obs'] = [['I3MCWeightDict.InteractionType',
                                       'NuMuProperties.muon_energy',
                                       'NuMuProperties.muon_energy_losses'],
                                      'Additional observables.']
    comp_dict_jm['n_files_total'] = [10000, 'Number of files in Madison.']
    comp_dict_jm['n_files_loaded'] = [[1000], 'Number of files loaded.']

    add_dict_to_traj(traj, 'components.mu_he', 'Mu Component.',
                     comp_dict_jm, to_log=False)

    #####
    # mu very high energy
    #####
    
    comp_dict_jmv = {}
    comp_dict_jmv['type'] = ['juliet_mu', 'Component type.']

    comp_dict_jmv['file_list'] = [
        [os.path.join(JULIET_PATH, 'mu_very_high_energy_merged_1kfiles.hdf5')],
        'Files for the component.']
    comp_dict_jmv['flux_names'] = [
        ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
         'ahlers_2010_1E18'],
        'Weight cols.']
    comp_dict_jmv['additional_obs'] = [['I3MCWeightDict.InteractionType',
                                       'NuMuProperties.muon_energy',
                                       'NuMuProperties.muon_energy_losses'],
                                      'Additional observables.']
    comp_dict_jmv['n_files_total'] = [2000, 'Number of files in Madison.']
    comp_dict_jmv['n_files_loaded'] = [[1000], 'Number of files loaded.']

    add_dict_to_traj(traj, 'components.mu_vhe', 'VHE Mu Component.',
                     comp_dict_jmv, to_log=False)

    #####
    # tau high energy
    #####
    
    comp_dict_jt = {}
    comp_dict_jt['type'] = ['juliet_tau', 'Component type.']

    comp_dict_jt['file_list'] = [
        [os.path.join(JULIET_PATH, 'tau_high_energy_merged_1kfiles.hdf5')],
        'Files for the component.']
    comp_dict_jt['flux_names'] = [
        ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
         'ahlers_2010_1E18'],
        'Weight cols.']
    comp_dict_jt['additional_obs'] = [['I3MCWeightDict.InteractionType',
                                       'NuTauProperties.tau_energy',
                                       'NuTauProperties.tau_energy_losses'],
                                      'Additional observables.']
    comp_dict_jt['n_files_total'] = [10000, 'Number of files in Madison.']
    comp_dict_jt['n_files_loaded'] = [[1000], 'Number of files loaded.']

    add_dict_to_traj(traj, 'components.tau_he', 'Tau Component.',
                     comp_dict_jt, to_log=False)

    #####
    # tau very high energy
    #####
    
    comp_dict_jtv = {}
    comp_dict_jtv['type'] = ['juliet_tau', 'Component type.']

    comp_dict_jtv['file_list'] = [
        [os.path.join(JULIET_PATH, 'tau_very_high_energy_merged_1kfiles.hdf5')],
        'Files for the component.']
    comp_dict_jtv['flux_names'] = [
        ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
         'ahlers_2010_1E18'],
        'Weight cols.']
    comp_dict_jtv['additional_obs'] = [['I3MCWeightDict.InteractionType',
                                       'NuTauProperties.muon_energy',
                                       'NuTauProperties.muon_energy_losses'],
                                      'Additional observables.']
    comp_dict_jtv['n_files_total'] = [2000, 'Number of files in Madison.']
    comp_dict_jtv['n_files_loaded'] = [[1000], 'Number of files loaded.']

    add_dict_to_traj(traj, 'components.tau_vhe', 'VHE Tau Component.',
                     comp_dict_jtv, to_log=False)

    #####
    # nutau high energy
    #####
    
    comp_dict_jnt = {}
    comp_dict_jnt['type'] = ['juliet_nutau', 'Component type.']

    comp_dict_jnt['file_list'] = [
        [os.path.join(JULIET_PATH, 'nutau_high_energy_merged_1kfiles.hdf5')],
        'Files for the component.']
    comp_dict_jnt['flux_names'] = [
        ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
         'ahlers_2010_1E18'],
        'Weight cols.']
    comp_dict_jnt['additional_obs'] = [['I3MCWeightDict.InteractionType',
                                        'NuTauProperties.tau_energy',
                                        'NuTauProperties.tau_energy_losses'],
                                      'Additional observables.']
    comp_dict_jnt['n_files_total'] = [10000, 'Number of files in Madison.']
    comp_dict_jnt['n_files_loaded'] = [[1000], 'Number of files loaded.']

    add_dict_to_traj(traj, 'components.nutau_he', 'NuTau Component.',
                     comp_dict_jnt, to_log=False)

    #####
    # nutau very high energy
    #####
    
    comp_dict_jntv = {}
    comp_dict_jntv['type'] = ['juliet_nutau', 'Component type.']

    comp_dict_jntv['file_list'] = [
        [os.path.join(JULIET_PATH, 'nutau_very_high_energy_merged_1kfiles.hdf5')],
        'Files for the component.']
    comp_dict_jntv['flux_names'] = [
        ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
         'ahlers_2010_1E18'],
        'Weight cols.']
    comp_dict_jntv['additional_obs'] = [['I3MCWeightDict.InteractionType',
                                       'NuTauProperties.tau_energy',
                                       'NuTauProperties.tau_energy_losses'],
                                      'Additional observables.']
    comp_dict_jntv['n_files_total'] = [2000, 'Number of files in Madison.']
    comp_dict_jntv['n_files_loaded'] = [[1000], 'Number of files loaded.']

    add_dict_to_traj(traj, 'components.nutau_vhe', 'VHE NuTau Component.',
                     comp_dict_jntv, to_log=False)

    #####
    # nue high energy
    #####
    
    comp_dict_jne = {}
    comp_dict_jne['type'] = ['juliet_nue', 'Component type.']

    comp_dict_jne['file_list'] = [
        [os.path.join(JULIET_PATH, 'nue_high_energy_merged_1kfiles.hdf5')],
        'Files for the component.']
    comp_dict_jne['flux_names'] = [
        ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
         'ahlers_2010_1E18'],
        'Weight cols.']
    comp_dict_jne['additional_obs'] = [['I3MCWeightDict.InteractionType',
                                        'NuEProperties.visible_energy'],
                                      'Additional observables.']
    comp_dict_jne['n_files_total'] = [10000, 'Number of files in Madison.']
    comp_dict_jne['n_files_loaded'] = [[1000], 'Number of files loaded.']

    add_dict_to_traj(traj, 'components.nue_he', 'NuE Component.',
                     comp_dict_jne, to_log=False)

    #####
    # nue very high energy
    #####
    
    comp_dict_jnev = {}
    comp_dict_jnev['type'] = ['juliet_nue', 'Component type.']

    comp_dict_jnev['file_list'] = [
        [os.path.join(JULIET_PATH, 'nue_very_high_energy_merged_1kfiles.hdf5')],
        'Files for the component.']
    comp_dict_jnev['flux_names'] = [
        ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
         'ahlers_2010_1E18'],
        'Weight cols.']
    comp_dict_jnev['additional_obs'] = [['I3MCWeightDict.InteractionType',
                                         'NuEProperties.visible_energy'],
                                      'Additional observables.']
    comp_dict_jnev['n_files_total'] = [2000, 'Number of files in Madison.']
    comp_dict_jnev['n_files_loaded'] = [[1000], 'Number of files loaded.']

    add_dict_to_traj(traj, 'components.nue_vhe', 'VHE NuE Component.',
                     comp_dict_jnev, to_log=False)


    
    return traj

def build_anatoli_trajcetory():

    traj = Trajectory(name='Trajectory', add_time=False)
    if os.uname().nodename in ['Brians-MBP', 'dhcp-172-16-174-147.private.pa.umd.edu']:
        DATA_PATH='/Users/brianclark/Documents/work/IceCube/ehe/new_i3files'
        JULIET_PATH = '/Users/brianclark/Documents/work/IceCube/ehe/new_i3files'
    elif os.uname().nodename in ['condor00']:
        DATA_PATH='/data/i3home/baclark/transfer_files'
        JULIET_PATH = '/data/i3home/baclark/transfer_files'
    else:
        raise NotImplementedError

    #####
    # numu high energy
    #####
    
    comp_dict_jnm = {}
    comp_dict_jnm['type'] = ['juliet_numu', 'Component type.']

    comp_dict_jnm['file_list'] = [
        [os.path.join(JULIET_PATH, 'numu_high_energy_merged_1kfiles.hdf5')],
        'Files for the component.']
    comp_dict_jnm['flux_names'] = [
        ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
         'ahlers_2010_1E18'],
        'Weight cols.']
    comp_dict_jnm['additional_obs'] = [['I3MCWeightDict.InteractionType',
                                       'NuMuProperties.muon_energy',
                                       'NuMuProperties.muon_energy_losses'],
                                      'Additional observables.']
    comp_dict_jnm['n_files_total'] = [10000, 'Number of files in Madison.']
    comp_dict_jnm['n_files_loaded'] = [[1000], 'Number of files loaded.']

    add_dict_to_traj(traj, 'components.numu_he', 'NuMu Component.',
                     comp_dict_jnm, to_log=False)


    # #####
    # # numu high energy
    # #####

    # comp_dict_jnmv = {}
    # comp_dict_jnmv['type'] = ['juliet_numu', 'Component type.']

    # comp_dict_jnmv['file_list'] = [
    #     [os.path.join(JULIET_PATH, 'numu_very_high_energy_merged_1kfiles.hdf5')],
    #     'Files for the component.']
    # comp_dict_jnmv['flux_names'] = [
    #     ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
    #      'ahlers_2010_1E18'],
    #     'Weight cols.']
    # comp_dict_jnmv['additional_obs'] = [['I3MCWeightDict.InteractionType',
    #                                    'NuMuProperties.muon_energy',
    #                                    'NuMuProperties.muon_energy_losses'],
    #                                   'Additional observables.']
    # comp_dict_jnmv['n_files_total'] = [2000, 'Number of files in Madison.']
    # comp_dict_jnmv['n_files_loaded'] = [[1000], 'Number of files loaded.']

    # add_dict_to_traj(traj, 'components.numu_vhe', 'VHE NuMu Component.',
    #                  comp_dict_jnmv, to_log=False)
    
    # #####
    # # mu high energy
    # #####
    
    # comp_dict_jm = {}
    # comp_dict_jm['type'] = ['juliet_mu', 'Component type.']

    # comp_dict_jm['file_list'] = [
    #     [os.path.join(JULIET_PATH, 'mu_high_energy_merged_1kfiles.hdf5')],
    #     'Files for the component.']
    # comp_dict_jm['flux_names'] = [
    #     ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
    #      'ahlers_2010_1E18'],
    #     'Weight cols.']
    # comp_dict_jm['additional_obs'] = [['I3MCWeightDict.InteractionType',
    #                                    'NuMuProperties.muon_energy',
    #                                    'NuMuProperties.muon_energy_losses'],
    #                                   'Additional observables.']
    # comp_dict_jm['n_files_total'] = [10000, 'Number of files in Madison.']
    # comp_dict_jm['n_files_loaded'] = [[1000], 'Number of files loaded.']

    # add_dict_to_traj(traj, 'components.mu_he', 'Mu Component.',
    #                  comp_dict_jm, to_log=False)

    # #####
    # # mu very high energy
    # #####
    
    # comp_dict_jmv = {}
    # comp_dict_jmv['type'] = ['juliet_mu', 'Component type.']

    # comp_dict_jmv['file_list'] = [
    #     [os.path.join(JULIET_PATH, 'mu_very_high_energy_merged_1kfiles.hdf5')],
    #     'Files for the component.']
    # comp_dict_jmv['flux_names'] = [
    #     ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
    #      'ahlers_2010_1E18'],
    #     'Weight cols.']
    # comp_dict_jmv['additional_obs'] = [['I3MCWeightDict.InteractionType',
    #                                    'NuMuProperties.muon_energy',
    #                                    'NuMuProperties.muon_energy_losses'],
    #                                   'Additional observables.']
    # comp_dict_jmv['n_files_total'] = [2000, 'Number of files in Madison.']
    # comp_dict_jmv['n_files_loaded'] = [[1000], 'Number of files loaded.']

    # add_dict_to_traj(traj, 'components.mu_vhe', 'VHE Mu Component.',
    #                  comp_dict_jmv, to_log=False)

    # #####
    # # tau high energy
    # #####
    
    # comp_dict_jt = {}
    # comp_dict_jt['type'] = ['juliet_tau', 'Component type.']

    # comp_dict_jt['file_list'] = [
    #     [os.path.join(JULIET_PATH, 'tau_high_energy_merged_1kfiles.hdf5')],
    #     'Files for the component.']
    # comp_dict_jt['flux_names'] = [
    #     ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
    #      'ahlers_2010_1E18'],
    #     'Weight cols.']
    # comp_dict_jt['additional_obs'] = [['I3MCWeightDict.InteractionType',
    #                                    'NuTauProperties.tau_energy',
    #                                    'NuTauProperties.tau_energy_losses'],
    #                                   'Additional observables.']
    # comp_dict_jt['n_files_total'] = [10000, 'Number of files in Madison.']
    # comp_dict_jt['n_files_loaded'] = [[1000], 'Number of files loaded.']

    # add_dict_to_traj(traj, 'components.tau_he', 'Tau Component.',
    #                  comp_dict_jt, to_log=False)

    # #####
    # # tau very high energy
    # #####
    
    # comp_dict_jtv = {}
    # comp_dict_jtv['type'] = ['juliet_tau', 'Component type.']

    # comp_dict_jtv['file_list'] = [
    #     [os.path.join(JULIET_PATH, 'tau_very_high_energy_merged_1kfiles.hdf5')],
    #     'Files for the component.']
    # comp_dict_jtv['flux_names'] = [
    #     ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
    #      'ahlers_2010_1E18'],
    #     'Weight cols.']
    # comp_dict_jtv['additional_obs'] = [['I3MCWeightDict.InteractionType',
    #                                    'NuTauProperties.muon_energy',
    #                                    'NuTauProperties.muon_energy_losses'],
    #                                   'Additional observables.']
    # comp_dict_jtv['n_files_total'] = [2000, 'Number of files in Madison.']
    # comp_dict_jtv['n_files_loaded'] = [[1000], 'Number of files loaded.']

    # add_dict_to_traj(traj, 'components.tau_vhe', 'VHE Tau Component.',
    #                  comp_dict_jtv, to_log=False)

    # #####
    # # nutau high energy
    # #####
    
    # comp_dict_jnt = {}
    # comp_dict_jnt['type'] = ['juliet_nutau', 'Component type.']

    # comp_dict_jnt['file_list'] = [
    #     [os.path.join(JULIET_PATH, 'nutau_high_energy_merged_1kfiles.hdf5')],
    #     'Files for the component.']
    # comp_dict_jnt['flux_names'] = [
    #     ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
    #      'ahlers_2010_1E18'],
    #     'Weight cols.']
    # comp_dict_jnt['additional_obs'] = [['I3MCWeightDict.InteractionType',
    #                                     'NuTauProperties.tau_energy',
    #                                     'NuTauProperties.tau_energy_losses'],
    #                                   'Additional observables.']
    # comp_dict_jnt['n_files_total'] = [10000, 'Number of files in Madison.']
    # comp_dict_jnt['n_files_loaded'] = [[1000], 'Number of files loaded.']

    # add_dict_to_traj(traj, 'components.nutau_he', 'NuTau Component.',
    #                  comp_dict_jnt, to_log=False)

    # #####
    # # nutau very high energy
    # #####
    
    # comp_dict_jntv = {}
    # comp_dict_jntv['type'] = ['juliet_nutau', 'Component type.']

    # comp_dict_jntv['file_list'] = [
    #     [os.path.join(JULIET_PATH, 'nutau_very_high_energy_merged_1kfiles.hdf5')],
    #     'Files for the component.']
    # comp_dict_jntv['flux_names'] = [
    #     ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
    #      'ahlers_2010_1E18'],
    #     'Weight cols.']
    # comp_dict_jntv['additional_obs'] = [['I3MCWeightDict.InteractionType',
    #                                    'NuTauProperties.tau_energy',
    #                                    'NuTauProperties.tau_energy_losses'],
    #                                   'Additional observables.']
    # comp_dict_jntv['n_files_total'] = [2000, 'Number of files in Madison.']
    # comp_dict_jntv['n_files_loaded'] = [[1000], 'Number of files loaded.']

    # add_dict_to_traj(traj, 'components.nutau_vhe', 'VHE NuTau Component.',
    #                  comp_dict_jntv, to_log=False)

    # #####
    # # nue high energy
    # #####
    
    # comp_dict_jne = {}
    # comp_dict_jne['type'] = ['juliet_nue', 'Component type.']

    # comp_dict_jne['file_list'] = [
    #     [os.path.join(JULIET_PATH, 'nue_high_energy_merged_1kfiles.hdf5')],
    #     'Files for the component.']
    # comp_dict_jne['flux_names'] = [
    #     ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
    #      'ahlers_2010_1E18'],
    #     'Weight cols.']
    # comp_dict_jne['additional_obs'] = [['I3MCWeightDict.InteractionType',
    #                                     'NuEProperties.visible_energy'],
    #                                   'Additional observables.']
    # comp_dict_jne['n_files_total'] = [10000, 'Number of files in Madison.']
    # comp_dict_jne['n_files_loaded'] = [[1000], 'Number of files loaded.']

    # add_dict_to_traj(traj, 'components.nue_he', 'NuE Component.',
    #                  comp_dict_jne, to_log=False)

    # #####
    # # nue very high energy
    # #####
    
    # comp_dict_jnev = {}
    # comp_dict_jnev['type'] = ['juliet_nue', 'Component type.']

    # comp_dict_jnev['file_list'] = [
    #     [os.path.join(JULIET_PATH, 'nue_very_high_energy_merged_1kfiles.hdf5')],
    #     'Files for the component.']
    # comp_dict_jnev['flux_names'] = [
    #     ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
    #      'ahlers_2010_1E18'],
    #     'Weight cols.']
    # comp_dict_jnev['additional_obs'] = [['I3MCWeightDict.InteractionType',
    #                                      'NuEProperties.visible_energy'],
    #                                   'Additional observables.']
    # comp_dict_jnev['n_files_total'] = [2000, 'Number of files in Madison.']
    # comp_dict_jnev['n_files_loaded'] = [[1000], 'Number of files loaded.']

    # add_dict_to_traj(traj, 'components.nue_vhe', 'VHE NuE Component.',
    #                  comp_dict_jnev, to_log=False)


    #####
    # corsika
    #####

    comp_dict_mu = {}
    comp_dict_mu['type'] = ['corsika', 'Component type.']
    comp_dict_mu['file_list'] = [
        [os.path.join(DATA_PATH, '20848_10k_files_merged.hd5'),
         os.path.join(DATA_PATH, '20848_10k_more_files_merged.hd5'),
         os.path.join(DATA_PATH, '22023_with_corrections/22023_0_1000.hd5'),
         os.path.join(DATA_PATH, '22023_with_corrections/22023_1_999.hd5'),
         os.path.join(DATA_PATH, '22023_with_corrections/22023_2_998.hd5'),
         os.path.join(DATA_PATH, '22023_with_corrections/22023_3_998.hd5'),
         os.path.join(DATA_PATH, '22023_with_corrections/22023_4_998.hd5'),
         os.path.join(DATA_PATH, '22023_with_corrections/22023_5_1000.hd5'),
         os.path.join(DATA_PATH, '22023_with_corrections/22023_6_1000.hd5'),
         os.path.join(DATA_PATH, '22023_with_corrections/22023_7_1000.hd5'),
         os.path.join(DATA_PATH, '22023_with_corrections/22023_8_1000.hd5'),
         os.path.join(DATA_PATH, '22023_with_corrections/22023_9_999.hd5'),
         os.path.join(DATA_PATH, '21962_with_corrections/21962_0_998.hd5'),
         os.path.join(DATA_PATH, '21962_with_corrections/21962_1_1000.hd5'),
         os.path.join(DATA_PATH, '21962_with_corrections/21962_2_997.hd5'),
         os.path.join(DATA_PATH, '21962_with_corrections/21962_3_995.hd5'),
         os.path.join(DATA_PATH, '21962_with_corrections/21962_4_999.hd5'),
         os.path.join(DATA_PATH, '21962_with_corrections/21962_5_999.hd5'),
         os.path.join(DATA_PATH, '21962_with_corrections/21962_6_997.hd5'),
         os.path.join(DATA_PATH, '21962_with_corrections/21962_7_996.hd5'),
         os.path.join(DATA_PATH, '21962_with_corrections/21962_8_992.hd5'),
         os.path.join(DATA_PATH, '21962_with_corrections/21962_9_996.hd5'),
         os.path.join(DATA_PATH, '22187_with_corrections/22187_0_982.hd5'),
         os.path.join(DATA_PATH, '22187_with_corrections/22187_1_989.hd5'),
         os.path.join(DATA_PATH, '22187_with_corrections/22187_2_987.hd5'),
         os.path.join(DATA_PATH, '22187_with_corrections/22187_3_986.hd5'),
         os.path.join(DATA_PATH, '22187_with_corrections/22187_4_987.hd5'),
         os.path.join(DATA_PATH, '22187_with_corrections/22187_5_986.hd5'),
         os.path.join(DATA_PATH, '22187_with_corrections/22187_6_981.hd5'),
         os.path.join(DATA_PATH, '22187_with_corrections/22187_7_982.hd5'),
         os.path.join(DATA_PATH, '22187_with_corrections/22187_8_982.hd5'),
         os.path.join(DATA_PATH, '22187_with_corrections/22187_9_979.hd5'),
        ],
        'Files for the component.']
    comp_dict_mu['flux_names'] = [
        ['GaisserH4a'],
        'Weight cols.']
    comp_dict_mu['n_files_total'] = [99812, 'Number of files in Madison.']
    comp_dict_mu['n_files_loaded'] = [
        [
            9981.2, 
            9981.2,
            
            1000., 
            999., 
            998., 
            998., 998., 1000., 1000., 1000., 1000., 999.,
            
            998., 
            1000., 
            997., 
            995., 999., 999., 997., 996., 992., 996.,
            
            982., 
            989., 
            987., 
            986., 987., 986., 981., 982., 982., 979.
        ],
        'Number of files loaded.']

    add_dict_to_traj(traj, 'components.corsika', 'Corsika Component.',
                     comp_dict_mu, to_log=False)
    
    return traj
    
def build_muon_test_trajectory():
    
    traj = Trajectory(name='Trajectory', add_time=False)
    if os.uname().nodename in ['Brians-MBP', 'dhcp-172-16-174-138.private.pa.umd.edu']:
        DATA_PATH='/Users/brianclark/Documents/work/IceCube/ehe/new_i3files'
        JULIET_PATH = '/Users/brianclark/Documents/work/IceCube/ehe/new_i3files'
    elif os.uname().nodename in ['condor00']:
        DATA_PATH='/data/i3home/baclark/transfer_files'
        JULIET_PATH = '/data/i3home/baclark/transfer_files'
    else:
        raise NotImplementedError
    
    j_flux_names = ['northern_tracks', 'northern_tracks_updated', 'ahlers_gzk',
         'ahlers_2010_1E18', 'cascades', 'hese_7']


    #####
    # numu high energy
    #####
    
    comp_dict_jnm = {}
    comp_dict_jnm['type'] = ['juliet_numu', 'Component type.']

    comp_dict_jnm['file_list'] = [
        [os.path.join(DATA_PATH, 'numu_high_energy_merged_1kfiles.hdf5')],
        'Files for the component.']
    comp_dict_jnm['flux_names'] = [
        j_flux_names,
        'Weight cols.']
    comp_dict_jnm['additional_obs'] = [['I3MCWeightDict.InteractionType',
                                       'NuMuProperties.muon_energy',
                                       'NuMuProperties.muon_energy_losses',
                                       'EHEMuMillipede_SplineMPEseed_depE.value',
                                       'EHEMuMillipede_SplineMPEseed_100m_depE.value',
                                       'EHEMuMillipede_SplineMPEseed_200m_depE.value',
                                       'EHEMuMillipede_SplineMPEseed_500m_depE.value', 
                                       ],
                                      'Additional observables.']
    comp_dict_jnm['n_files_total'] = [2000, 'Number of files in Madison.']
    comp_dict_jnm['n_files_loaded'] = [[1000], 'Number of files loaded.']

    add_dict_to_traj(traj, 'components.numu_he', 'NuMu Component.',
                     comp_dict_jnm, to_log=False)

    #####
    # numu very high energy
    #####

    comp_dict_jnmv = {}
    comp_dict_jnmv['type'] = ['juliet_numu', 'Component type.']

    comp_dict_jnmv['file_list'] = [
        [os.path.join(DATA_PATH, 'numu_very_high_energy_merged_1kfiles.hdf5')],
        'Files for the component.']
    comp_dict_jnmv['flux_names'] = [
        j_flux_names,
        'Weight cols.']
    comp_dict_jnmv['additional_obs'] = [['I3MCWeightDict.InteractionType',
                                       'NuMuProperties.muon_energy',
                                       'NuMuProperties.muon_energy_losses',
                                       'EHEMuMillipede_SplineMPEseed_depE.value',
                                       'EHEMuMillipede_SplineMPEseed_100m_depE.value',
                                       'EHEMuMillipede_SplineMPEseed_200m_depE.value',
                                       'EHEMuMillipede_SplineMPEseed_500m_depE.value',
                                        ],
                                      'Additional observables.']
    comp_dict_jnmv['n_files_total'] = [2000, 'Number of files in Madison.']
    comp_dict_jnmv['n_files_loaded'] = [[1000], 'Number of files loaded.']

    add_dict_to_traj(traj, 'components.numu_vhe', 'VHE NuMu Component.',
                     comp_dict_jnmv, to_log=False)
        
    comp_dict_mu = {}
    comp_dict_mu['type'] = ['corsika', 'Component type.']
    comp_dict_mu['file_list'] = [
        [os.path.join(DATA_PATH, '20848_10k_files_merged.hd5'),
        #  os.path.join(DATA_PATH, '20848_10k_more_files_merged.hd5'),
         os.path.join(DATA_PATH, '22023_with_corrections/22023_0_1000.hd5'),
        #  os.path.join(DATA_PATH, '22023_with_corrections/22023_1_999.hd5'),
        #  os.path.join(DATA_PATH, '22023_with_corrections/22023_2_998.hd5'),
        #  os.path.join(DATA_PATH, '22023_with_corrections/22023_3_998.hd5'),
        #  os.path.join(DATA_PATH, '22023_with_corrections/22023_4_998.hd5'),
        #  os.path.join(DATA_PATH, '22023_with_corrections/22023_5_1000.hd5'),
        #  os.path.join(DATA_PATH, '22023_with_corrections/22023_6_1000.hd5'),
        #  os.path.join(DATA_PATH, '22023_with_corrections/22023_7_1000.hd5'),
        #  os.path.join(DATA_PATH, '22023_with_corrections/22023_8_1000.hd5'),
        #  os.path.join(DATA_PATH, '22023_with_corrections/22023_9_999.hd5'),
         os.path.join(DATA_PATH, '21962_with_corrections/21962_0_998.hd5'),
        #  os.path.join(DATA_PATH, '21962_with_corrections/21962_1_1000.hd5'),
        #  os.path.join(DATA_PATH, '21962_with_corrections/21962_2_997.hd5'),
        #  os.path.join(DATA_PATH, '21962_with_corrections/21962_3_995.hd5'),
        #  os.path.join(DATA_PATH, '21962_with_corrections/21962_4_999.hd5'),
        #  os.path.join(DATA_PATH, '21962_with_corrections/21962_5_999.hd5'),
        #  os.path.join(DATA_PATH, '21962_with_corrections/21962_6_997.hd5'),
        #  os.path.join(DATA_PATH, '21962_with_corrections/21962_7_996.hd5'),
        #  os.path.join(DATA_PATH, '21962_with_corrections/21962_8_992.hd5'),
        #  os.path.join(DATA_PATH, '21962_with_corrections/21962_9_996.hd5'),
         os.path.join(DATA_PATH, '22187_with_corrections/22187_0_982.hd5'),
        #  os.path.join(DATA_PATH, '22187_with_corrections/22187_1_989.hd5'),
        #  os.path.join(DATA_PATH, '22187_with_corrections/22187_2_987.hd5'),
        #  os.path.join(DATA_PATH, '22187_with_corrections/22187_3_986.hd5'),
        #  os.path.join(DATA_PATH, '22187_with_corrections/22187_4_987.hd5'),
        #  os.path.join(DATA_PATH, '22187_with_corrections/22187_5_986.hd5'),
        #  os.path.join(DATA_PATH, '22187_with_corrections/22187_6_981.hd5'),
        #  os.path.join(DATA_PATH, '22187_with_corrections/22187_7_982.hd5'),
        #  os.path.join(DATA_PATH, '22187_with_corrections/22187_8_982.hd5'),
        #  os.path.join(DATA_PATH, '22187_with_corrections/22187_9_979.hd5'),
        ],
        'Files for the component.']
    comp_dict_mu['flux_names'] = [
        ['GaisserH4a'],
        'Weight cols.']
    comp_dict_mu['n_files_total'] = [99812, 'Number of files in Madison.']
    comp_dict_mu['n_files_loaded'] = [
        [
            9981.2, 
            # 9981.2,
            1000., 
            # 999., 998., 998., 998., 1000., 1000., 1000., 1000., 999.,
            998., 
            # 1000., 997., 995., 999., 999., 997., 996., 992., 996.,
            982., 
            # 989., 987., 986., 987., 986., 981., 982., 982., 979.
        ],
        'Number of files loaded.']

    add_dict_to_traj(traj, 'components.corsika', 'Corsika Component.',
                     comp_dict_mu, to_log=False)
    
    return traj

def build_angres_trajectory():
    traj = Trajectory(name='Trajectory', add_time=False)
    if os.uname().nodename in ['Brians-MBP', 'dhcp-172-16-174-138.private.pa.umd.edu']:
        DATA_PATH='/Users/brianclark/Documents/work/IceCube/ehe/new_i3files'
        JULIET_PATH = '/Users/brianclark/Documents/work/IceCube/ehe/new_i3files'
    elif os.uname().nodename in ['condor00']:
        DATA_PATH='/data/i3home/baclark/transfer_files'
        JULIET_PATH = '/data/i3home/baclark/transfer_files'
    else:
        raise NotImplementedError
    
    # numu events from nugen (cuz why not?)
    comp_dict_nm = {}
    comp_dict_nm['type'] = ['numu', 'Component type.']

    comp_dict_nm['file_list'] = [
        [os.path.join(DATA_PATH, '21220_merged.hd5')],
        'Files for the component.']
    comp_dict_nm['flux_names'] = [
        ['northern_tracks'],
        'Weight cols.']
    comp_dict_nm['additional_obs'] = [['I3MCWeightDict.InteractionType',
                                       'NuMuProperties.muon_energy',
                                       'NuMuProperties.muon_energy_losses'],
                                      'Additional observables.']
    comp_dict_nm['n_files_total'] = [10000, 'Number of files in Madison.']
    comp_dict_nm['n_files_loaded'] = [[10000], 'Number of files loaded.']

    add_dict_to_traj(traj, 'components.numu', 'NuMu Component.',
                     comp_dict_nm, to_log=False)

    
    # corsika
    comp_dict_mu = {}
    comp_dict_mu['type'] = ['corsika', 'Component type.']
    comp_dict_mu['file_list'] = [
        [os.path.join(DATA_PATH, '20848_10k_files_merged.hd5'),
        #  os.path.join(DATA_PATH, '20848_10k_more_files_merged.hd5'),
         os.path.join(DATA_PATH, '22023_with_corrections/22023_0_1000.hd5'),
        #  os.path.join(DATA_PATH, '22023_with_corrections/22023_1_999.hd5'),
        #  os.path.join(DATA_PATH, '22023_with_corrections/22023_2_998.hd5'),
        #  os.path.join(DATA_PATH, '22023_with_corrections/22023_3_998.hd5'),
        #  os.path.join(DATA_PATH, '22023_with_corrections/22023_4_998.hd5'),
        #  os.path.join(DATA_PATH, '22023_with_corrections/22023_5_1000.hd5'),
        #  os.path.join(DATA_PATH, '22023_with_corrections/22023_6_1000.hd5'),
        #  os.path.join(DATA_PATH, '22023_with_corrections/22023_7_1000.hd5'),
        #  os.path.join(DATA_PATH, '22023_with_corrections/22023_8_1000.hd5'),
        #  os.path.join(DATA_PATH, '22023_with_corrections/22023_9_999.hd5'),
         os.path.join(DATA_PATH, '21962_with_corrections/21962_0_998.hd5'),
        #  os.path.join(DATA_PATH, '21962_with_corrections/21962_1_1000.hd5'),
        #  os.path.join(DATA_PATH, '21962_with_corrections/21962_2_997.hd5'),
        #  os.path.join(DATA_PATH, '21962_with_corrections/21962_3_995.hd5'),
        #  os.path.join(DATA_PATH, '21962_with_corrections/21962_4_999.hd5'),
        #  os.path.join(DATA_PATH, '21962_with_corrections/21962_5_999.hd5'),
        #  os.path.join(DATA_PATH, '21962_with_corrections/21962_6_997.hd5'),
        #  os.path.join(DATA_PATH, '21962_with_corrections/21962_7_996.hd5'),
        #  os.path.join(DATA_PATH, '21962_with_corrections/21962_8_992.hd5'),
        #  os.path.join(DATA_PATH, '21962_with_corrections/21962_9_996.hd5'),
         os.path.join(DATA_PATH, '22187_with_corrections/22187_0_982.hd5'),
        #  os.path.join(DATA_PATH, '22187_with_corrections/22187_1_989.hd5'),
        #  os.path.join(DATA_PATH, '22187_with_corrections/22187_2_987.hd5'),
        #  os.path.join(DATA_PATH, '22187_with_corrections/22187_3_986.hd5'),
        #  os.path.join(DATA_PATH, '22187_with_corrections/22187_4_987.hd5'),
        #  os.path.join(DATA_PATH, '22187_with_corrections/22187_5_986.hd5'),
        #  os.path.join(DATA_PATH, '22187_with_corrections/22187_6_981.hd5'),
        #  os.path.join(DATA_PATH, '22187_with_corrections/22187_7_982.hd5'),
        #  os.path.join(DATA_PATH, '22187_with_corrections/22187_8_982.hd5'),
        #  os.path.join(DATA_PATH, '22187_with_corrections/22187_9_979.hd5'),
        ],
        'Files for the component.']
    comp_dict_mu['flux_names'] = [
        ['GaisserH4a'],
        'Weight cols.']
    comp_dict_mu['n_files_total'] = [99812, 'Number of files in Madison.']
    comp_dict_mu['n_files_loaded'] = [
        [
            9981.2, 
            # 9981.2,
            1000., 
            # 999., 998., 998., 998., 1000., 1000., 1000., 1000., 999.,
            998., 
            # 1000., 997., 995., 999., 999., 997., 996., 992., 996.,
            982., 
            # 989., 987., 986., 987., 986., 981., 982., 982., 979.
        ],
        'Number of files loaded.']

    add_dict_to_traj(traj, 'components.corsika', 'Corsika Component.',
                     comp_dict_mu, to_log=False)
    
    return traj

to_ignore=[
        'Borderness',
        'Dustyness',
        'EHEMillipedeQuantiles',
        'EHEMillipedeQuantiles_MCPrimary',
        'EHEMuMillipede_MCPrimaryseed_vecd',
        'EHE_SplineMPECharacteristics',
        'EHE_SplineMPEDirectHitsC',
        'EHE_SplineMPEFitParams',
        'EHE_SplineMPE_TruncatedEnergy_SpiceMie_AllBINS_Muon', 
        'EHE_SplineMPE_TruncatedEnergy_SpiceMie_AllBINS_Neutrino',
        'EHE_SplineMPE_TruncatedEnergy_SpiceMie_AllDOMS_Muon',
        'EHE_SplineMPE_TruncatedEnergy_SpiceMie_AllDOMS_Neutrino',
        'EHE_SplineMPE_TruncatedEnergy_SpiceMie_BINS_Muon',
        'EHE_SplineMPE_TruncatedEnergy_SpiceMie_BINS_Neutrino',
        'EHE_SplineMPE_TruncatedEnergy_SpiceMie_BINS_dEdxVector',
        'EHE_SplineMPE_TruncatedEnergy_SpiceMie_ORIG_Muon',
        'EHE_SplineMPE_TruncatedEnergy_SpiceMie_ORIG_Neutrino',
        'EnergyLossLabels_vecd',
        'EnergyLossLabels_Truth_vecd',
        'LabelsQuantiles',
        'LabelsQuantiles_Truth',
        'MPEFitParaboloid',
        'MPEFitParaboloidFitParams',
        'MuProperties',
        'TauProperties',
        'NuMuProperties',
        'NuEProperties',
        'NuTauProperties',
        'PolyplopiaPrimary',
        'ProjectedQ',
        'TruncatedStochasticity',
        'ddddr_100CascadeParams',
        'ddddr_100Params',
        'ddddr_100_features',
        'ddddr_150CascadeParams',
        'ddddr_150Params',
        'ddddr_150_features'
]
        

def get_aeff(resource_name):
    '''
    Return energy bin centers and effective area
    
    Energy bin centers should be in units of GeV
    Effective area should be in units of meters^2 * steradian
    
    '''
    
    energy_GeV = None
    aeffsr_m2sr = None
    
    if resource_name == 'ehenextgen':
        data = np.load('ehenextgen_total_aeff_scale.npz')
        energy_GeV = data['bin_center']
        aeffsr_m2sr =data['avg_aeff_m2sr']
    else:
        raise NotImplementedError

    return energy_GeV, aeffsr_m2sr

def get_aeff_vs_zen(resource_name):
    
    energy_bins_GeV = None
    czen_bins = None
    aeff_m2 = None
    
    if resource_name == 'ehenextgen':
        data = np.load('ehenextgen_total_aeff_vs_e_and_zen_scale.npz')
        energy_bins_GeV = data['energy_bins']
        czen_bins = data['czen_bins']
        aeff_m2 = data['avg_aeff_m2']
    else:
        raise NotImplementedError

    return energy_bins_GeV, czen_bins, aeff_m2
    
    
def find_nearest_energy_bin(array, value):
    array = np.asarray(array)
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]

def get_bin_centers(bins):
    return (bins[1:] + bins[:-1]) * 0.5

def get_bin_widths(bins):
     return bins[:-1] - bins[1:]
