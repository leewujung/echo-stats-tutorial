
#############################################################################
## 2016/08/10
### Mixed assemblage phasor summation
The following codes were revised from the old codes on the WHOI computer for new calculation:

- `simulation_mixed_assem_loop_20160809`

- `fig_mixed_assem_ka_cmp`

<figure>
   <img src=".\figs_results\fig_mixed_assem_ka_cmp\fig_mixed_assem_ka_cmp_Ns01_smpl1e7.png" width="600">
</figure>

- `fig_mixed_assem_Ns_cmp`

<figure>
   <img src=".\figs_results\fig_mixed_assem_Ns_cmp\fig_mixed_assem_Ns_cmp_ka14pi_smpl1e7.png" width="600">
</figure>

- `fig_mixed_assem_set_cmp`

Seems like 1e7 is necessary to obtain the "convergence" while 1e5 or 1e6 samples are not enough.

<figure>
   <img src=".\figs_results\fig_mixed_assem_set_cmp\fig_mixed_assem_set_cmp_Ns01_ka14pi_smplN1e7.png" width="600">
</figure>
<figure>
   <img src=".\figs_results\fig_mixed_assem_set_cmp\fig_mixed_assem_set_cmp_Ns01_ka14pi_smplN1e6.png" width="600">
</figure>
<figure>
   <img src=".\figs_results\fig_mixed_assem_set_cmp\fig_mixed_assem_set_cmp_Ns01_ka14pi_smplN1e5.png" width="600">
</figure>


- `fig_mixed_assem_smplN_cmp`

<figure>
   <img src=".\figs_results\fig_mixed_assem_smplN_cmp\fig_mixed_assem_smplN_cmp_Ns01_ka14pi.png" width="600">
</figure>
<figure>
   <img src=".\figs_results\fig_mixed_assem_smplN_cmp\fig_mixed_assem_smplN_cmp_Ns10_ka14pi.png" width="600">
</figure>
<figure>
   <img src=".\figs_results\fig_mixed_assem_smplN_cmp\fig_mixed_assem_smplN_cmp_Ns50_ka14pi.png" width="600">
</figure>



#############################################################################
## 2016/10/19
### Plot Fig. 5 and Fig. 12 for tutorial
Codes:

- `fig_5_rayleigh_linlog` -- PDF, CDF, and PFA for Rayleigh distribution
- `fig_12_pb_ka` -- Beampattern PDF (Pb) for beamwidths (1,3,5,10 degs)

Spent some time on the routine to find the good-enough ka values for approaching the exact 2-way *full* beamwidth. It is in the first part of code `fig_12_pb_ka`.

