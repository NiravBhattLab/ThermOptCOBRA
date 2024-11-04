# ThermOptCOBRA
Codes for "Thermodynamically Consistent Construction and Analysis of Genome-Scale Metabolic Models"   

Authors: Pavan Kumar S and Nirav P Bhatt
<p align="center">
  <img src="https://github.com/NiravBhattLab/ThermOptiCOBRA/blob/main/PaperThermOptCOBRA/Figures/BioRenderFigures/OverviewOfTOCS.png" alt="ThermOptiCS" width="500"/>
</p>

**PaperThermOptCOBRA:**
Codes required to recreate all the results and figures shown in the article

**ThermOptCOBRA:**
Matlab codes of the algorithms developed in this work

# ThermOptEnumMILP: To enumerate all the TICs in a metabolic model
<dl class="mat function">
<dt class="sig sig-object mat" id="ThermOptEnumMILP">
<span class="sig-name descname"><span class="pre">ThermOptEnumMILP</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="pre">model</span></em>, <em class="sig-param"><span class="pre">timeLimit</span></em><span class="sig-paren">)</span></a></dt>
<dd><p>Enumerates all the Thermodynamically infeasible cycles in a given model</p>
<dl class="field-list simple">
<dt class="field-odd">USAGE<span class="colon">:</span></dt>
<dd class="field-odd"><p><strong>[TICs,Direction,TIC_Rxns,modModel,opt] = ThermOptEnumMILP</strong> (<em>model,timeLimit</em>)</p>
</dd>
<dt class="field-even">INPUTS<span class="colon">:</span></dt>
<dd class="field-even"><p><strong>model</strong> – COBRA model structure for which TICs has be found</p>
</dd>
<dt class="field-odd">OPTIONAL INPUTS<span class="colon">:</span></dt>
<dd class="field-odd"><p><strong>timeLimit</strong> – If the algorithm takes more than this time null set will
be returned for all the outputs</p>
</dd>
<dt class="field-even">OUTPUTS<span class="colon">:</span></dt>
<dd class="field-even"><ul class="simple">
<li><p><strong>TICs</strong> – List of all the Thermodynamically infeasible cycles in
the given input model</p></li>
<li><p><strong>Direction</strong> – Relative flux coefficients for reactions
in the corresponding TICs</p></li>
<li><p><strong>TIC_Rxns</strong> – Reaction list that participates in the TICs</p></li>
<li><p><strong>modModel</strong> – Modified model that has no irreversible reactions that
carry flux in reverse direction</p></li>
<li><p><strong>opt</strong> – Says whether the provided solution is optimal or not</p></li>
</ul>
</dd>
</dl>
</dd></dl>

# ThermOptCC: To identify thermodynamically feasible flux directions

<dl class="mat function">
<dt class="sig sig-object mat" id="ThermOptCC">
<span class="sig-name descname"><span class="pre">ThermOptCC</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="pre">model</span></em>, <em class="sig-param"><span class="pre">tol</span></em>, <em class="sig-param"><span class="pre">TICs</span></em>, <em class="sig-param"><span class="pre">Dir</span></em><span class="sig-paren">)</span></a></dt>
<dd><p>Identifies thermodynamically feasible flux directions for all the
reactions in the input model</p>
<dl class="field-list simple">
<dt class="field-odd">USAGE<span class="colon">:</span></dt>
<dd class="field-odd"><p><strong>[a,TICs,Dir] = ThermOptCC</strong> (<em>model,tol</em>)</p>
</dd>
<dt class="field-even">INPUTS<span class="colon">:</span></dt>
<dd class="field-even"><ul class="simple">
<li><p><strong>model</strong> – COBRA model structure for which thermodynamic feasibility
of the reactions has to be identified</p></li>
<li><p><strong>tol</strong> – Tolerance value (User defined non-zero value).</p></li>
</ul>
</dd>
<dt class="field-odd">OPTIONAL INPUTS<span class="colon">:</span></dt>
<dd class="field-odd"><ul class="simple">
<li><p><strong>TICs</strong> – List of all the Thermodynamically infeasible cycles in
the given input model</p></li>
<li><p><strong>Dir</strong> – The flux directions for reactions in the corresponding
TICs</p></li>
</ul>
</dd>
<dt class="field-even">OUTPUTS<span class="colon">:</span></dt>
<dd class="field-even"><ul class="simple">
<li><p><strong>a</strong> – A cell describing the thermodynamically feasible direction
of the reactions in the given input model</p></li>
<li><p><strong>TICs</strong> – List of all the Thermodynamically infeasible cycles in
the given input model</p></li>
<li><p><strong>Dir</strong> – The flux directions for reactions in the corresponding
TICs</p></li>
</ul>
</dd>
</dl>
</dd></dl>

# ThermOptiCS: To build thermodynamically consistent CSM

<dl class="mat function">
<dt class="sig sig-object mat" id="ThermOptiCS">
<span class="sig-name descname"><span class="pre">ThermOptiCS</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="pre">model</span></em>, <em class="sig-param"><span class="pre">core</span></em>, <em class="sig-param"><span class="pre">tol</span></em>, <em class="sig-param"><span class="pre">TICs</span></em>, <em class="sig-param"><span class="pre">Dir</span></em><span class="sig-paren">)</span></a></dt>
<dd><p>Thermodynamically feasible context-specific model building.</p>
<dl class="field-list simple">
<dt class="field-odd">USAGE<span class="colon">:</span></dt>
<dd class="field-odd"><p><strong>[Model,bCoreRxns,TICs,Dir] = ThermOptiCS</strong> (<em>model,core,tol</em>)</p>
</dd>
<dt class="field-even">INPUTS<span class="colon">:</span></dt>
<dd class="field-even"><ul class="simple">
<li><p><strong>model</strong> – COBRA model structure for which act as the GSMM</p></li>
<li><p><strong>core</strong> – Reaction IDs which are defined to be core (These
reactions will be present in the final model)</p></li>
<li><p><strong>tol</strong> – Tolerance value (User defined non-zero value).</p></li>
</ul>
</dd>
<dt class="field-odd">OUTPUTS<span class="colon">:</span></dt>
<dd class="field-odd"><ul class="simple">
<li><p><strong>Model</strong> – Context-specific model</p></li>
<li><p><strong>bCoreRxns</strong> – Core reactions that are thermodynamically blocked (or infeasible)</p></li>
<li><p><strong>TICs</strong> – List of all the Thermodynamically infeasible cycles in
the given input model</p></li>
<li><p><strong>Dir</strong> – The flux directions for reactions in the corresponding
TICs</p></li>
</ul>
</dd>
</dl>
</dd></dl>

# ThermOptFlux: To obtain TIC-free flux

<dl class="mat function">
<dt class="sig sig-object mat" id="ThermOptFlux">
<span class="sig-name descname"><span class="pre">ThermOptFlux</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="pre">model</span></em>, <em class="sig-param"><span class="pre">flux</span></em>, <em class="sig-param"><span class="pre">TICs</span></em>, <em class="sig-param"><span class="pre">Dir</span></em><span class="sig-paren">)</span></a></dt>
<dd><p>Gets the TIC free flux values for the given model and a feasible flux</p>
<dl class="field-list simple">
<dt class="field-odd">USAGE<span class="colon">:</span></dt>
<dd class="field-odd"><p><strong>[flux] = ThermOptFlux</strong> (<em>model,flux,TICs,Dir</em>)</p>
</dd>
<dt class="field-even">INPUTS<span class="colon">:</span></dt>
<dd class="field-even"><ul class="simple">
<li><p><strong>model</strong> – COBRA model structure for from which the flux values are
obtained</p></li>
<li><p><strong>flux</strong> – Flux distribution obtained from any flux analysis methods
(FBA or flux sampling)</p></li>
<li><p><strong>TICs</strong> – List of all the Thermodynamically infeasible cycles in
the given input model</p></li>
<li><p><strong>Dir</strong> – The flux directions for reactions in the corresponding
TICs</p></li>
</ul>
</dd>
<dt class="field-odd">OUTPUTS<span class="colon">:</span></dt>
<dd class="field-odd"><p><strong>flux</strong> – TIC free flux</p>
</dd>
</dl>
</dd></dl>


# ThemOptEnumLP (This code is not tested. Recommended to use ThermOptEnumMILP.m)

<dl class="mat function">
<dt class="sig sig-object mat" id="ThermOptEnumLP">
<span class="sig-name descname"><span class="pre">ThermOptEnumLP</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="pre">model</span></em><span class="sig-paren">)</span></a></dt>
<dd><p>Enumerates all the Thermodynamically infeasible cycles in a given model</p>
<dl class="field-list simple">
<dt class="field-odd">USAGE<span class="colon">:</span></dt>
<dd class="field-odd"><p><strong>[TICs,Direction,TIC_Rxns,modModel] = ThermOptEnumLP</strong> (<em>model</em>)</p>
</dd>
<dt class="field-even">INPUTS<span class="colon">:</span></dt>
<dd class="field-even"><p><strong>model</strong> – COBRA model structure for which TICs has be found</p>
</dd>
<dt class="field-odd">OUTPUTS<span class="colon">:</span></dt>
<dd class="field-odd"><ul class="simple">
<li><p><strong>TICs</strong> – List of all the Thermodynamically infeasible cycles in
the given input model</p></li>
<li><p><strong>Direction</strong> – Relative flux coefficients for reactions
in the corresponding TICs</p></li>
<li><p><strong>TIC_Rxns</strong> – Reaction list that participates in the TICs</p></li>
<li><p><strong>modModel</strong> – Modified model that has no irreversible reactions that
carry flux in reverse direction</p></li>
</ul>
</dd>
</dl>
</dd></dl>

