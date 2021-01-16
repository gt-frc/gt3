f#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 14:58:56 2018

@author: max
"""
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('font', **{'family':'sans-serif', 'sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font', **{'family':'serif', 'serif':['Palatino']})
rc('text', usetex=True)
import sys

class gt3plots():
    sys.dont_write_bytecode = True 
    #if ITER model
    #if myshot.param.d3d_iter==1:
    #    fig = plt.figure(figsize=(8, 8))
    #    ax1 = fig.add_subplot(111)
    #    ax1.set_title(r'ITER First Wall', fontsize=20, )
    #    ax1.set_ylabel(r'Z ($m$)', fontsize=20)
    #    ax1.set_xlabel(r'R ($m$)', fontsize=20)
    #    ax1.axis('equal')
    #    ax1.plot(myshot.exp.lim_ITER[:, 0], myshot.exp.lim_ITER[:, 1])
    
    #print ('test = ', myshot.ntrl.sol_flx_line[0])
    
    mill_plot = plt.figure(figsize=(12, 8))
    ax1 = mill_plot.add_subplot(111)
    ax1.axis('equal')
    ax1.set_xlim(3.5, 7)
    ax1.set_ylim(-4.5, -2.5)
    ax1.contour(myshot.brnd.R, myshot.brnd.Z, myshot.brnd.r, 30)
    ax1.plot(myshot.brnd.R[-1, :], myshot.brnd.Z[-1, :])
    sep = np.vstack((myshot.exp.bdry, myshot.exp.bdry[0]))
    #ax1.plot(myshot.ntrl.points_unique[:, 0], myshot.ntrl.points_unique[:, 1], 'o')
    #ax1.plot(myshot.ntrl.sol_pts[:, 0], myshot.ntrl.sol_pts[:, 1], color='black')
    for i in range(0, 5):
        #ax1.plot(myshot.ntrl.sol_flx_line[:, 2*i], myshot.ntrl.sol_flx_line[:, 2*i+1], color='black')
        ax1.plot(myshot.ntrl.sol_pts[i*100:(i+1)*100, 0], myshot.ntrl.sol_pts[i*100:(i+1)*100, 1], color='black')
    ax1.plot(myshot.exp.lim_vertex[:, 0], myshot.exp.lim_vertex[:, 1])
    #ax1.plot(myshot.ntrl.sol_lim_pts[:, 0], myshot.ntrl.sol_lim_pts[:, 1], 'o', color='red')
    #ax1.plot(sep[:, 0], sep[:, 1], color='green', lw=3)
    #CS = ax1.contourf(myshot.brnd.R, myshot.brnd.Z, myshot.brnd.B_p, 500) #plot something calculated by miller
    #plt.colorbar(CS)
    
    
    ITER_wall_x = np.array([4.04688, 4.04688, 4.30729, 4.98438, 5.73958, 7.38021, 7.88802, 
                            8.23958, 8.39583, 8.30469, 7.90104, 7.28906, 6.27344, 6.33854, 
                            6.27344, 6.15625, 6.05208, 5.96094, 5.86979, 5.77865, 5.71354, 
                            5.66146, 5.60938, 5.57031, 5.55729, 5.55729, 5.27083, 5.25781, 
                            5.25781, 5.20573, 5.15365, 5.07552, 4.95833, 4.85417, 4.6849, 
                            4.64583, 4.5026, 4.1901, 4.35938, 4.48958, 4.51562, 4.48958, 
                            4.42448, 4.33333, 4.21615, 4.125, 4.04688])
    ITER_wall_y = np.array([-2.56193, 3.60495, 4.28292, 4.66102, 4.53064, 3.16167, 2.43155, 
                            1.66232, 0.59322, -0.436767, -1.34941, -2.2751, -3.05737, 
                            -3.16167, -3.22686, -3.22686, -3.2399, -3.27901, -3.33116, 
                            -3.39635, -3.46154, -3.53977, -3.63103, -3.83963, -3.95698, 
                            -4.53064, -4.26988, -4.16558, -4.00913, -3.91786, -3.85267, 
                            -3.78748, -3.73533, -3.72229, -3.73533, -3.76141, -3.91786, 
                            -3.90482, -3.52673, -3.21382, -3.08344, -2.9661, -2.78357, 
                            -2.67927, -2.61408, -2.57497, -2.56193])
    
    R0 = myshot.param.R0_a
    a = myshot.param.a
    mill_plot2 = plt.figure(figsize=(8, 12))
    ax1 = mill_plot2.add_subplot(111)
    #ax1.set_title(r'Example Miller Flux Surfaces', fontsize=30)
    #ax1.set_ylabel(r'Z ($m$)', fontsize=30)
    #ax1.set_xlabel(r'R ($m$)', fontsize=30)
    #ax1.tick_params(labelsize=30)
    ax1.axis('equal')
    #print myshot.brnd.R.T[::-2]
    for i, vals in enumerate(myshot.brnd.R.T):
        if i%4 == 0:
            ax1.plot(myshot.brnd.R.T[i], myshot.brnd.Z.T[i], 'gray', alpha=0.5, lw=1)
    #ax1.plot(myshot.brnd.R[-1], myshot.brnd.Z[-1], 'black')
    CS1 = ax1.contour(myshot.brnd.R, myshot.brnd.Z, myshot.brnd.rho, levels = [0.333, 0.666, 0.999], colors='black', linewidths=6)
    
    fmt = {}
    strs = [r'$\rho = 0.33$', r'$\rho = 0.66$', r'$\rho = 1$']
    for l, s in zip(CS1.levels, strs):
        fmt[l] = s
    
    manual_locations = [(R0+a/3, 0.1), (R0+2*a/3, -0.2), (R0+a, -0.5)]
    # Label every other level using strings
    #plt.clabel(CS1, CS1.levels[::], inline=True, fmt=fmt, fontsize=10)
    
    labels = ax1.clabel(CS1, CS1.levels[::], inline=1, fmt=fmt, fontsize=20, manual = manual_locations)
    for l in labels:
        l.set_rotation(0)
    print('myshot.param.d3d_iter = ', myshot.param.d3d_iter)
    if myshot.param.d3d_iter==1:
        ax1.plot(ITER_wall_x, ITER_wall_y, lw=6)
    else:
        ax1.plot(myshot.exp.lim_vertex[:, 0], myshot.exp.lim_vertex[:, 1], lw=6)
    #ax1.plot(myshot.exp.lim_ITER[:, 0], myshot.exp.lim[:, 1], 'red')
    ax1.set_axis_off()
    mill_plot2.tight_layout()
    mill_plot2.savefig("/home/max/Nextcloud/Max/max_phd/Codes/X-MILLER/example_miller.png", dpi=300, transparent=True) 
    
    
    fig = plt.figure(figsize=(12, 8))
    fig.suptitle(r'Using original stuff')
    ax1 = fig.add_subplot(231)
    ax1.set_title(r'ni', fontsize=20)
    ax1.set_ylabel(r'ni', fontsize=20)
    ax1.set_xlabel(r'rho', fontsize=20)
    ax1.plot(myshot.brnd.rho[:, 0], myshot.brnd.ni[:, 0])
    
    ax2 = fig.add_subplot(232)
    ax2.set_title(r'ne', fontsize=20)
    ax2.set_ylabel(r'ne', fontsize=20)
    ax2.set_xlabel(r'rho', fontsize=20)
    ax2.plot(myshot.brnd.rho[:, 0], myshot.brnd.ne[:, 0])
    
    ax3 = fig.add_subplot(233)
    ax3.set_title(r'Ti', fontsize=20)
    ax3.set_ylabel(r'Ti', fontsize=20)
    ax3.set_xlabel(r'rho', fontsize=20)
    ax3.plot(myshot.brnd.rho[:, 0], myshot.brnd.Ti_kev[:, 0])
    
    ax4 = fig.add_subplot(234)
    ax4.set_title(r'Te', fontsize=20, )
    ax4.set_ylabel(r'Te', fontsize=20)
    ax4.set_xlabel(r'rho', fontsize=20)
    ax4.plot(myshot.brnd.rho[:, 0], myshot.brnd.Te_kev[:, 0])
    plt.tight_layout()
    
    ax4 = fig.add_subplot(235)
    ax4.set_title(r'E_pot', fontsize=20, )
    ax4.set_ylabel(r'E_pot', fontsize=20)
    ax4.set_xlabel(r'rho', fontsize=20)
    ax4.plot(myshot.brnd.rho[:, 0], myshot.brnd.E_pot[:, 0])
    plt.tight_layout()
    
    ax4 = fig.add_subplot(236)
    ax4.set_title(r'B_phi', fontsize=20, )
    ax4.set_ylabel(r'B_phi', fontsize=20)
    ax4.set_xlabel(r'R', fontsize=20)
    ax4.plot(myshot.brnd.R[:, myshot.brnd.R.shape[1]/2], myshot.brnd.B_phi[:, myshot.brnd.R.shape[1]/2], 'blue')
    ax4.plot(myshot.brnd.R[:, 0], myshot.brnd.B_phi[:, 0], 'blue')
    plt.tight_layout()
    
    #ax1.plot(myshot.brnd.R[-1, :], myshot.brnd.Z[-1, :])
    #sep = np.vstack((myshot.exp.bdry, myshot.exp.bdry[0]))
    #ax1.plot(sep[:, 0], sep[:, 1], color='green', lw=3)
    #CS = ax1.pcolor(myshot.brnd.R, myshot.brnd.Z, myshot.brnd.B_p, vmin=0, vmax=0.6) #plot something calculated by miller
    #CS = ax1.pcolor(myshot.brnd.R, myshot.brnd.Z, myshot.brnd.Psi_norm) #plot something calculated by miller
    #CS = ax1.pcolor(myshot.brnd.R, myshot.brnd.Z, B_p_in, )
    #plt.colorbar(CS)
    """
    """
    B_p_in = griddata(np.column_stack((myshot.exp.R.flatten(), myshot.exp.Z.flatten())), 
                            myshot.exp.B_p.flatten(), 
                            (myshot.brnd.R, myshot.brnd.Z), 
                            method='cubic')
    
    
    mill_plot3 = plt.figure(figsize=(8, 8))
    ax1 = mill_plot3.add_subplot(111)
    ax1.axis('equal')
    #ax1.contour(myshot.brnd.R, myshot.brnd.Z, myshot.brnd.r, 70)
    #ax1.plot(myshot.brnd.R[-1, :], myshot.brnd.Z[-1, :])
    #sep = np.vstack((myshot.exp.bdry, myshot.exp.bdry[0]))
    #ax1.plot(sep[:, 0], sep[:, 1], color='green', lw=3)
    CS = ax1.contourf(myshot.brnd.R, myshot.brnd.Z, myshot.brnd.B_p, 500)
    #CS = ax1.pcolor(myshot.brnd.R, myshot.brnd.Z, B_p_in, vmin=0, vmax=0.6)
    #CS = ax1.contourf(myshot.brnd.R, myshot.brnd.Z, B_p_in, 500, vmin=0, vmax=0.6) #plot something calculated by miller
    plt.colorbar(CS)
    
    
    
    """
    fig = plt.figure(figsize=(6, 6))
    #if myshot.param.d3d_iter==1:
        #fig.suptitle('Calculated IOL in ITER', fontsize=24, y=0.8)
    #else:
        #fig.suptitle('Calculated IOL in DIII-D', fontsize=24, y=0.8)
    ax1 = fig.add_subplot(111)
    #ax1.set_title(r'Calculated $E_{orb}$ in DIII-D', fontsize=20)
    ax1.set_title(r'Increase $E_{orb}$', fontsize=20)
    #ax1.set_ylabel(r'Cumulative energy loss fraction', fontsize=20)
    ax1.set_ylabel(r'Increase in cumulative energy loss', fontsize=20)
    ax1.set_xlabel(r'$\rho$', fontsize=20)
    #ax1.text(0.05, 0.85, 
    #    "Baseline Rad. = 53.492 MW\n$T_i$ increase Rad. = %s MW"%(round(rad_tot_adpak/1E6, 3)), 
    #    fontsize=15, 
    #    transform = ax1.transAxes, 
    #    bbox={'facecolor':'white', 'alpha':1, 'pad':10})
    
    ax1.grid(b=True, which='both', axis='both')
    ax1.set_xlim(0.75, 1)
    #ax1.set_ylim(0, 1)
    #ax1.plot(np.linspace(0, 1, myshot.param.rpts), myshot.tiol.E_orb, label='Eorb', color='black', lw=3)
    #ax1.plot(np.linspace(0, 1, myshot.param.rpts), E_orb_ITER, label='Baseline ITER', color='black', lw=3)
    #ax1.plot(np.linspace(0, 1, myshot.param.rpts), E_orb_ITER_T2, '--', label='20\% Temp increase', color='red', lw=3)
    
    ax1.plot(np.linspace(0, 1, myshot.param.rpts), (E_orb_ITER_T2 - E_orb_ITER), label='Increase', color='black', lw=3)
    
    ax1.legend(bbox_to_anchor=(0.015, 0.6), 
        loc='lower left', 
        ncol=1, 
        prop={'size':20}, 
        #transform = ax1.transAxes, 
        shadow=True)  
    """
    
    #fig2 = plt.figure(figsize=(6, 6))
    #if myshot.param.d3d_iter==1:
        #fig.suptitle('Calculated IOL in ITER', fontsize=24, y=0.8)
    #else:
        #fig.suptitle('Calculated IOL in DIII-D', fontsize=24, y=0.8)
    #ax1 = fig2.add_subplot(111)
    #ax1.contourf(myshot.brnd.R[:450, :], myshot.brnd.Z[:450, :], myshot.brnd.rho[:450, :], 1, colors='blue')
    #ax1.contourf(myshot.brnd.R[450:, :], myshot.brnd.Z[450:, :], myshot.brnd.rho[450:, :], 1, colors='green')
    #ax1.plot(myshot.brnd.R[-1, :], myshot.brnd.Z[-1, :], 'black', lw=3)
    #ax1.axis('equal')
    #plt.axis('off')
    
    
    """
    ax2 = fig.add_subplot(132)
    ax2.set_title(r'$M_{orb}$', fontsize=20)
    ax2.set_ylabel(r'cumulative momentum loss fraction', fontsize=20)
    ax2.set_xlabel(r'minor radius', fontsize=20)
    ax2.grid(b=True, which='both', axis='both')
    ax2.set_xlim(0.75, 1)
    ax2.plot(np.linspace(0, 1, myshot.param.rpts), myshot.tiol.M_orb, label='Morb', color='black', lw=3)
    
    ax3 = fig.add_subplot(133)
    ax3.set_title(r'$E_{orb}$', fontsize=20)
    ax3.set_ylabel(r'cumulative energy loss fraction', fontsize=20)
    ax3.set_xlabel(r'minor radius', fontsize=20)
    ax3.grid(b=True, which='both', axis='both')
    ax3.set_xlim(0.75, 1)
    ax3.plot(np.linspace(0, 1, myshot.param.rpts), myshot.tiol.E_orb, label='Eorb', color='black', lw=3)
    plt.tight_layout()
    fig.subplots_adjust(top=0.84)
    """
    """
    myshot = plasma('toplasma1')
    myshot.solve()
    np.savetxt('input1_fiol_t.dat', myshot.tiol.F_orb, delimiter=', ')
    np.savetxt('input1_miol_t.dat', myshot.tiol.M_orb, delimiter=', ')
    np.savetxt('input1_eiol_t.dat', myshot.tiol.E_orb, delimiter=', ')
    
    myshot = plasma('toplasma2')
    myshot.solve()
    np.savetxt('input2_fiol_t.dat', myshot.tiol.F_orb, delimiter=', ')
    np.savetxt('input2_miol_t.dat', myshot.tiol.M_orb, delimiter=', ')
    np.savetxt('input2_eiol_t.dat', myshot.tiol.E_orb, delimiter=', ')
    
    myshot = plasma('toplasma3')
    myshot.solve()
    np.savetxt('input3_fiol_t.dat', myshot.tiol.F_orb, delimiter=', ')
    np.savetxt('input3_miol_t.dat', myshot.tiol.M_orb, delimiter=', ')
    np.savetxt('input3_eiol_t.dat', myshot.tiol.E_orb, delimiter=', ')
    """