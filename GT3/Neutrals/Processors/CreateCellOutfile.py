#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pandas as pd


def create_cell_outfile(neutpy_inst, outfile, midpts):
    df = pd.DataFrame()
    df['R'] = pd.Series(midpts[:,0], name='R')
    df['Z'] = pd.Series(midpts[:,1], name='Z')
    df['n_n_slow'] = pd.Series(neutpy_inst.nn.s, name='n_n_slow')
    df['n_n_thermal'] = pd.Series(neutpy_inst.nn.t, name='n_n_thermal')
    df['n_n_total'] = pd.Series(neutpy_inst.nn.tot, name='n_n_total')
    df['izn_rate_slow'] = pd.Series(neutpy_inst.izn_rate.s, name='izn_rate_slow')
    df['izn_rate_thermal'] = pd.Series(neutpy_inst.izn_rate.t, name='izn_rate_thermal')
    df['izn_rate_total'] = pd.Series(neutpy_inst.izn_rate.tot, name='izn_rate_total')
    #cell_df = iterate_namedtuple(neut.cell, df)
    df.to_csv(outfile)
