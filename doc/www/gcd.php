<a name="gcd"></a> 
<div class="divName"> 
<h2 class="name">Name:</h2> <?php linkTo("command","gcd","gcd");?> 
<span class="smallDescription">Computes the greatest common divider of polynomials or numbers. 
</span> 
</div> 
<div class="divLibraryName"> 
<h2 class="libraryname">Library name:</h2> 
<span class="commandline type">sollya_obj_t sollya_lib_gcd(sollya_obj_t, sollya_obj_t)</span> 
</div> 
<div class="divUsage"> 
<h2 class="category">Usage: </h2> 
<span class="commandline"><?php linkTo("command","gcd","gcd");?>(<span class="arg">p</span>, <span class="arg">q</span>) : (<span class="type">function</span>, <span class="type">function</span>) -&gt; <span class="type">function</span></span> 
 
</div> 
<div class="divParameters"> 
<h2 class="category">Parameters: </h2> 
<ul> 
<li><span class="arg">p</span> is a polynomial.</li> 
<li><span class="arg">q</span> is a polynomial.</li> 
</ul> 
</div> 
<div class="divDescription"> 
<h2 class="category">Description: </h2><ul> 
<li>TODO 
</ul> 
</div> 
<div class="divExamples"> 
<div class="divExample"> 
<h2 class="category">Example 1: </h2> 
&nbsp;&nbsp;&nbsp;&gt; gcd(x^2 + 2 * x + 1, x + 1);<br> 
&nbsp;&nbsp;&nbsp;1 + x<br> 
</div> 
</div> 
<div class="divSeeAlso"> 
<span class="category">See also: </span><?php linkTo("command","eucldiv","div");?>, <?php linkTo("command","euclmod","mod");?>, <?php linkTo("command","numberroots","numberroots");?> 
</div> 
