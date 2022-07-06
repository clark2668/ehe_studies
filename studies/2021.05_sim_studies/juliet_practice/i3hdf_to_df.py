#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division, absolute_import
from __future__ import with_statement, print_function
import glob
import warnings
import numpy as np
import tables
import pandas as pd
try:
    from tqdm import tqdm
except ImportError:
    NO_TQDM = True
else:
    NO_TQDM = False


class ObservableName(object):
    def __init__(self,
                 table_name=None,
                 col_name=None,
                 obs_name=None):

        assert ((table_name is not None and col_name is not None) or
                (obs_name is not None)), \
               ('Either \'table_name\' + \'col_name\' or' +
                ' \'obs_name\' must be provided')
        if table_name is None or col_name is None:
            self.tab, self.col = self.__split_obs_str(obs_name)
        else:
            self.tab = table_name
            self.col = col_name
        self.id_col_log = None

    def __str__(self):
        return '%s.%s' % (self.tab, self.col)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, othr):
        return ((self.tab, self.col) == (othr.tab, othr.col))

    def __ne__(self, othr):
        return ~self.__eq__(othr)

    def __hash__(self):
        return hash((self.tab, self.col))

    def __lt__(self, othr):
        return self.__str__() < othr.__str__()

    def __split_obs_str(self, obs):
        splitted = obs.split('.')
        assert len(splitted) >= 2, \
            'Falsy obbservable string syntax! Needs to be <tab>.<col>.'
        if len(splitted) > 2:
            col = splitted[-1]
            tab = str.join('.', splitted[:-1])
        elif len(splitted) == 2:
            col = splitted[1]
            tab = splitted[0]
        return tab, col


def split_obs_str(obs):
    if not isinstance(obs, list):
        obs = [obs]
    obs_dict = {}
    for o in obs:
        splitted = o.split('.')
        key = splitted[0]
        current_content = obs_dict.get(key, [])
        current_content.append(splitted[1])
        obs_dict[key] = current_content
    return obs_dict


def get_values_from_table(table, cols, dtype=float):
    values = np.empty((table.nrows, len(cols)), dtype=float)
    for i, row in enumerate(table.iterrows()):
        values[i, :] = [row[col] for col in cols]
    return values


class HDFcontainer:
    def __init__(self,
                 file_list=None,
                 directory=None,
                 id_cols=['Run',
                          'Event',
                          'SubEvent'],
                 exists_col=None,
                 silent=False):
        assert (directory is not None) or (file_list is not None), \
            'If component is not from aggregation directory or file_list'\
            'is needed!'
        if file_list is None:
            if '*' not in directory:
                directory += '*'
            file_list = glob.glob(directory)
        self.file_list = file_list
        self.id_cols = id_cols
        self.exists_col = exists_col
        self.id_col_log = None
        self.silent = silent

    def get_observables(self,
                        blacklist_tabs=[],
                        blacklist_cols=[],
                        blacklist_obs=[],
                        check_all=False,
                        return_str=False):
        component_set = set()
        if check_all:
            files = self.file_list
        else:
            files = self.file_list[:1]
        blacklist_cols.extend(self.id_cols)
        if self.exists_col is not None:
            blacklist_cols.append(self.exists_col)
        for i, file_name in enumerate(files):
            file_set = set()
            f = tables.open_file(file_name)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                for table in f.iter_nodes('/', classname='Table'):
                    table_name = table.name
                    if table_name not in blacklist_tabs:
                        col_names = [c for c in table.colnames
                                     if c not in blacklist_cols]
                        for c in col_names:
                            obs = ObservableName(table_name=table_name,
                                                 col_name=c)
                            if obs not in blacklist_obs:
                                file_set.add(obs)
                f.close()
                if i == 0:
                    component_set = file_set
                else:
                    component_set = component_set.intersection(file_set)
        return sorted(component_set)

    def get_values(self, table_key, cols):
        if isinstance(cols, str):
            cols = [cols]
        for i, file_name in enumerate(self.file_list):
            with warnings.catch_warnings():
                f = pd.HDFStore(file_name, 'r')
            table = f[table_key]
            drops = [c for c in table.columns
                     if c not in cols and c not in self.id_cols and
                     c != self.exists_col]
            table.drop(drops, axis=1, inplace=True)
            table.set_index(self.id_cols, inplace=True)
            f.close()
            id_log_series = pd.Series(data=i, index=table.index)
            if i == 0:
                if self.id_col_log is None:
                    id_col_log = id_log_series
                values = table
            else:
                if self.id_col_log is None:
                    id_col_log = id_col_log.append(id_log_series)
                values = values.append(table)
        if self.exists_col is not None:
            mask = values.get(self.exists_col) == 0
            values[mask] = np.NaN
            values.drop(self.exists_col, axis=1, inplace=True)
        rename_dict = {col: '%s.%s' % (table_key, col)
                       for col in values.columns}
        values.rename(columns=rename_dict, inplace=True)
        if self.id_col_log is None:
            self.id_col_log = id_col_log
        return values

    def create_obs_dict(self, observables):
        obs_dict = {}
        for obs in observables:
            if not isinstance(obs, ObservableName):
                try:
                    obs = ObservableName(obs_name=obs)
                except AttributeError:
                    pass
            if obs.tab in list(obs_dict.keys()):
                obs_dict[obs.tab].append(obs.col)
            else:
                obs_dict[obs.tab] = [obs.col]
        return obs_dict

    def get_df(self, observables):
        obs_dict = self.create_obs_dict(observables)
        n_obs = len(observables)
        if self.silent or NO_TQDM:
            for i, [table_key, cols] in enumerate(obs_dict.items()):
                tab_values = self.get_values(table_key, cols)
                if i == 0:
                    df = tab_values
                else:
                    df = df.join(tab_values, how='outer')
        else:
            with tqdm(total=n_obs, unit='Observables') as pbar:
                for i, [table_key, cols] in enumerate(obs_dict.items()):
                    try:
                        tab_values = self.get_values(table_key, cols)
                    except KeyError:
                        print('Skipping to load {}!'.format(table_key))
                        continue
                    for col in cols:
                        if col.endswith('vector_index'):
                            tab_values = self.unduplicate_table(tab_values,
                                                                table_key,
                                                                col)
                    if i == 0:
                        df = tab_values
                    else:
                        df = df.join(tab_values, how='outer')
                    pbar.update(len(cols))
        return df

    def unduplicate_table(self, df, tab_key, col):
        print(f'Unduplicating Table for {tab_key}.{col}')
        col = '.'.join([tab_key, col])
        max_index = int(np.max(df[col]))
        item_col = col.replace('vector_index', 'item')

        def multiindex_pivot(df, index=None, columns=None, values=None):
            if index is None:
                names = list(df.index.names)
                df = df.reset_index()
            else:
                names = index
            list_index = df[names].values
            tuples_index = [tuple(i) for i in list_index] # hashable
            df = df.assign(tuples_index=tuples_index)
            df = df.pivot(index="tuples_index",
                          columns=columns,
                          values=values
            )
            tuples_index = df.index  # reduced
            index = pd.MultiIndex.from_tuples(tuples_index,
                                              names=names)
            df.index = index
            return df

        unduplicated_df = df.pipe(
            multiindex_pivot,
            index=None,
            columns=col,
            values=item_col
        )

        new_colname_prefix = col.replace('vector_index', 'vector_elem')
        unduplicated_df.columns = [new_colname_prefix + f'_{int(i):d}'
                                   for i in unduplicated_df.columns.values]

        # for i in range(0, max_index+1):
        #     df[col.replace('vector_index', f'vector_elem_{i}')] = np.nan

        # mask_0 = df[col] == 0
        # for i in range(0, max_index+1):
        #     mask = df[col] == i
        #     df.loc[mask_0, col.replace('vector_index', f'vector_elem_{i}')] = \
        #         df.loc[mask, item_col]

        # unduplicated_df = df[~df.index.duplicated(keep='first')]
        return unduplicated_df

    def setup_table_header(self, df, needed_cols):
        dtypes = []
        names = []
        with tables.open_file(self.file_list[0], 'r') as f:
            node_iter = f.iter_nodes('/', classname='Table')
            for i in needed_cols:
                dtype = ref_node.dtype[i]
                dtypes.append((i, dtype))
                names.append(i)
        for i, dtype in zip(df.columns, df.dtypes):
            dtypes.append((str(i), dtype))
            names.append(i)
        return np.dtype(dtypes), names

    def append_table(self,
                     name,
                     df,
                     fill_dict={'SubEventStream': 0}):
        origninal_cols = df.columns
        needed_cols = [str(i) for i in self.id_cols]
        fill_cols = [str(i) for i in list(fill_dict.keys())]
        if self.exists_col is not None:
            needed_cols.append(str(self.exists_col))
        needed_cols.extend(fill_cols)
        dtypes, col_names = self.setup_table_header(df, needed_cols)
        for i, val in fill_dict.items():
            df[i] = val
        df[self.exists_col] = 1
        for i, file_name in enumerate(self.file_list):
            idx = self.id_col_log.loc[self.id_col_log == i].index
            temp_df = df.loc[idx]
            temp_df.fillna(0, inplace=True)
            temp_df = temp_df.reset_index()
            rec_arr = temp_df.loc[:, col_names].to_records(index=False)
            with tables.open_file(file_name, 'a') as f:
                new_tab = f.create_table('/', name, dtypes)
                correct_arr = np.rec.fromrecords(rec_arr, dtype=new_tab.dtype)
                new_tab.append(correct_arr)
        added_cols = [i for i in df.columns if i not in origninal_cols]
        df = df.drop(added_cols, axis=1)
