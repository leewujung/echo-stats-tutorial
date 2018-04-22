
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



**************************************************************
## 2016/10/19
### Plot Fig. 5 and Fig. 12 for tutorial
Codes:

- `fig_5_rayleigh_linlog` -- PDF, CDF, and PFA for Rayleigh distribution
- `fig_12_pb_ka` -- Beampattern PDF (Pb) for beamwidths (1,3,5,10 degs)

Spent some time on the routine to find the good-enough ka values for approaching the exact 2-way *full* beamwidth. It is in the first part of code `fig_12_pb_ka`.


**************************************************************
## 2017/04/17
### Prolate spheroid echo pdf (Fig. 17)
1. Update code to allow 3D rotation of prolate spheroid
2. Compare 2D and 3D results

### Update bbechopdf code
1. Use struct to pass input parameters... the code is MUCH cleaner now!
2. Still need work:
	- prolate spheroid option
	- get_tx function, now generate an octave bandwidth chirp centered at 50 kHz
	- restrained beampattern angle: variable 'bpa'


**************************************************************
## 2018/01/24
### Working on response to reviewers
- C.  Utility of the physics based models for inference: In Sect. II.A the authors make an important point that physics based models are useful for prediction (the forward problem) and for inference (the backward problem or estimation). I did not see the inference topic explored in the latter part of the text in any detail and fear this application may suffer from the same issues as matched field processing. The authors should at a minimum describe/reference how much of this has been done in terms of real applications and performance limits like Cramer Rao lower bounds.
  - Tim: I am not sure which section the reviewer was referring to. I couldn't find related materials in Sec. II.A. But I have put a few sentences below that could be morphed to integrate with the flow better somewhere... perhaps in Sec. III.A?? I only wrote about what we did, and didn't put in anything from the seafloor/sea surface domain. I also didn't mention what could be done -- that doesn't seem to belong here since it is oo detailed and not the focus of this tutorial.
  - **Draft response**: The majority of work reviewed in this tutorial focuses on the prediction of the echo PDF using a physics-based approach, and therefore is more in line with the "forward" problem. However, as explained above, the forward, predictive model is of great importance while performing inference for information of the scattering processes (the "inverse" problem). Specifically, the echo PDF model of multiple types and sizes of scatterers (Sec. VII.C.) has been applied to analyze simulated data (Lee and Stanton, 2014). The echo PDF model from pulsed broadband signals (Sec. VIII.C.) was applied to analyze fish echo data collected from the ocean (Lee and Stanton 2015). These studies demonstrated the significance of employing models with matching physical scattering processes in inference problems, and the capability of statistics-based analysis (which does not require absolute calibration) in providing inference results similar to those from energy-based models (which requires absolute calibration).


- 404: PFA carries with it the assumption of the noise-only hypothesis in a detection algorithm. What you’ve described here is the “exceedance distribution function” (=1-CDF). If you keep the PFA nomenclature, you need to describe it in the context of a detection test. I suggest calling it an EDF and noting it might be PFA if the echoes are unwanted or PD if the echoes are from the desired object.
  - **Previous email 1**: My understanding of EDF ([empirical distribution function]( https://en.wikipedia.org/wiki/Empirical_distribution_function)) is that it is used to refer to the distribution derived from an empirical sample. While the EDF can converge to the true distribution, they are not really equivalent. It is true though that by displaying PDF/PFAs from numerical simulations, we are using supposedly converged EDF to substitute for the true entity. But for nomenclature I don't agree with this substitution. I can see what the reviewer is trying to emphasize, but I think that's exactly why this quantity is called PFA, because we don't know a priori whether the signal we see that exceed a certain level is from a target or noise...
  - **Previous email 2**: We thought it's better to not use the [EDF]( https://en.wikipedia.org/wiki/Empirical_distribution_function) terminology since it is "empirical" and usually used to refer to the empirical estimate of CDF. As for PFA vs PD, we thought it would be the best to add a few sentences to explain the PFA vs PD context (as shown in [this figure](https://goo.gl/images/bbm8Po) when the PFA terminology is first introduced, but in the remain text keep the PFA terminology.
  - **Draft response 2018/01/25**: We thank the reviewer for the very informative discussion. We decided to keep using PFA in the manuscript, because in typical practical cases of using statistics of echoes to find target, the forward PDF models we introduced in this manuscript are useful in determining when the amplitude of an echo actually exceed an "expected" amplitude and should be classified as potential targets. In these cases we do not know a priori whether the signal is from a target, or we do not have enough information about the echo properties of desirable targets to be detected. Therefore using PD seems to be misleading. However we do appreciate the reviewer's suggestion and have added a sentence discussing how the quantity 1-CDF could be interpreted as PD (Probability of Detection) in some cases.


- 568: CFs can and are applied to multi-dimensional random variables
  - **Draft response**: This is true. We have removed the misleading sentence.


- 575-577: Suggest you mention the alternative of a numerical inversion of the CFs via the Hankel transform to directly obtain the envelope PDF, which avoids your infinite series summation. The Hankel-transform integral relationship between the CF of a circularly symmetric complex RV & the envelope PDF is in [D.M. Drumheller, “Pade approximations to matched filter amplitude probability functions ,” IEEE T-AES, 35:3, 1033-1045, 1999] & there are various numerical routines for Hankel transforms.
  - **Draft response**: We thank the reviewer for pointing out this source of reference. We have now added a sentence to mention the use of Hankel transform to compute the envelope PDF, with references.
  - Reference:
    - Drumheller, D. M. (1999). Pade approximations to matched filter amplitude probability functions. Aerospace and Electronic Systems, IEEE Transactions on, 35(3), 1033–1045. https://doi.org/10.1109/7.784072
    - Drumheller, D.M.; Lew, H. (2002). Homodyned-K fluctuation model. Aerospace and Electronic Systems, IEEE Transactions on, 38(2), 527–542. https://doi.org/10.1109/TAES.2002.1008984


- 1489: for the randomly oriented prolate spheroid, did you simulate the response using a randomly drawn orientation to dictate the corresponding power for a Rayleigh-distributed amplitude or did you do random draws from the PDF described by (28)?
  - **Draft response**: The simulation was done by randomly drawing an orientation of the prolate spheroid and use the output of eq.(24) multiplied with a Rayleigh-distributed modulation for simulating rough surfaces (as discussed in Sec. VI.B.3.c) to produce the echo amplitude of 1 sample. This process was repeated 1e7 times to arrive at an ensemble, based on which the various PDFs were estimated.


- 1497: important to be specific about how many trials you use in the Monte Carlo so the accuracy of the estimate can be assessed.
  - **Draft response**: We have added the number of samples in the Monte Carlo simulation in each figures.
  - Tim: I think we have discussed last time to add this.


- 1499 Fig. 15: the noise seen at low values of the PDFs here and in other figures should be explained.
  - **Draft response**: The noise seen at the low values of the PDFs are due to the choice of binning the data in log scale. As a results, within the same bin in log scale, the bins in linear scale at lower values are smaller than the bins at higher values. The "noise" seen on the estimated PDF results from the low number of Monte Carlo samples falling into the smaller bins (in linear scale) at low echo amplitude values.


- 1740 Sect. VIII.A: shouldn’t there be an “auto-correlation function PDF” as there is a “beampattern PDF”? If not, why not?
  - **Draft response**: We are not sure what the reviewer meant by an "auto-correlation function PDF". We assumed that the reviewer was referring to the matched-filter processing stage shown in Fig. 22. It is true that we can estimate a PDF based on the auto-correlated (perhaps more properly denoted as cross-correlation) function. However, because the incorporation of scatterer response and system response were conducted in the frequency domain and then transformed-sampled in the time domain, unless these two functions are both delta functions in time (uniform frequency response), the auto-correlation/cross-correlation function PDF does not yield representative echo PDFs that are more useful in real-world applications.
