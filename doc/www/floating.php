<a name="floating"></a> 
<div class="divName"> 
<h2 class="name">Name:</h2> <?php linkTo("command","floating","floating");?> 
<span class="smallDescription">indicates that floating-point formats should be used for <?php linkTo("command","fpminimax","fpminimax");?> 
</span> 
</div> 
<div class="divLibraryName"> 
<h2 class="libraryname">Library names:</h2> 
<span class="commandline type">sollya_obj_t sollya_lib_floating()</span> 
<span class="commandline type">int sollya_lib_is_floating(sollya_obj_t)</span> 
</div> 
<div class="divUsage"> 
<h2 class="category">Usage: </h2> 
<span class="commandline"><?php linkTo("command","floating","floating");?> : <span class="type">fixed|floating</span></span> 
 
</div> 
<div class="divDescription"> 
<h2 class="category">Description: </h2><ul> 
<li>The use of <?php linkTo("command","floating","floating");?> in the command <?php linkTo("command","fpminimax","fpminimax");?> indicates that the list of 
formats given as argument is to be considered to be a list of floating-point 
formats. 
See <?php linkTo("command","fpminimax","fpminimax");?> for details. 
</ul> 
</div> 
<div class="divExamples"> 
<div class="divExample"> 
<h2 class="category">Example 1: </h2> 
&nbsp;&nbsp;&nbsp;&gt; fpminimax(cos(x),6,[|D...|],[-1;1],floating);<br> 
&nbsp;&nbsp;&nbsp;0.99999974816012215939053930924274027347564697265625 + x * (-2.795931796958502334440230695107655659202089892465e-15 + x * (-0.49999286980201401719980935922649223357439041137695 + x * (4.0484539189054105169841244454207387920433372507922e-14 + x * (4.1633515528919168291466235132247675210237503051758e-2 + x * (-4.015858818743733758578949218474363725507386355118e-14 + x * (-1.3382240885483781024645200119493892998434603214264e-3))))))<br> 
</div> 
</div> 
<div class="divSeeAlso"> 
<span class="category">See also: </span><?php linkTo("command","fpminimax","fpminimax");?>, <?php linkTo("command","fixed","fixed");?> 
</div> 
