<a name="dirtysimplify"></a> 
<div class="divName"> 
<h2 class="name">Name:</h2> <?php linkTo("command","dirtysimplify","dirtysimplify");?> 
<span class="smallDescription">simplifies an expression representing a function 
</span> 
</div> 
<div class="divLibraryName"> 
<h2 class="libraryname">Library name:</h2> 
<span class="commandline type">sollya_obj_t sollya_lib_dirtysimplify(sollya_obj_t)</span> 
</div> 
<div class="divUsage"> 
<h2 class="category">Usage: </h2> 
<span class="commandline"><?php linkTo("command","dirtysimplify","dirtysimplify");?>(<span class="arg">function</span>) : <span class="type">function</span> -&gt; <span class="type">function</span></span> 
 
</div> 
<div class="divParameters"> 
<h2 class="category">Parameters: </h2> 
<ul> 
<li><span class="arg">function</span> represents the expression to be simplified</li> 
</ul> 
</div> 
<div class="divDescription"> 
<h2 class="category">Description: </h2><ul> 
<li>The command <?php linkTo("command","dirtysimplify","dirtysimplify");?> simplifies constant subexpressions of the 
expression given in argument representing the function 
<span class="arg">function</span>. Those constant subexpressions are evaluated using 
floating-point arithmetic with the global precision <?php linkTo("command","prec","prec");?>. 
</ul> 
</div> 
<div class="divExamples"> 
<div class="divExample"> 
<h2 class="category">Example 1: </h2> 
&nbsp;&nbsp;&nbsp;&gt; print(dirtysimplify(sin(pi * x)));<br> 
&nbsp;&nbsp;&nbsp;sin(3.1415926535897932384626433832795028841971693993751 * x)<br> 
&nbsp;&nbsp;&nbsp;&gt; print(dirtysimplify(erf(exp(3) + x * log(4))));<br> 
&nbsp;&nbsp;&nbsp;erf(20.0855369231876677409285296545817178969879078385544 + x * 1.3862943611198906188344642429163531361510002687205)<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 2: </h2> 
&nbsp;&nbsp;&nbsp;&gt; prec = 20!;<br> 
&nbsp;&nbsp;&nbsp;&gt; t = erf(0.5);<br> 
&nbsp;&nbsp;&nbsp;&gt; s = dirtysimplify(erf(0.5));<br> 
&nbsp;&nbsp;&nbsp;&gt; prec = 200!;<br> 
&nbsp;&nbsp;&nbsp;&gt; t;<br> 
&nbsp;&nbsp;&nbsp;0.520499877813046537682746653891964528736451575757963700058806<br> 
&nbsp;&nbsp;&nbsp;&gt; s;<br> 
&nbsp;&nbsp;&nbsp;0.52050018310546875<br> 
</div> 
</div> 
<div class="divSeeAlso"> 
<span class="category">See also: </span><?php linkTo("command","simplify","simplify");?>, <?php linkTo("command","autosimplify","autosimplify");?>, <?php linkTo("command","prec","prec");?>, <?php linkTo("command","evaluate","evaluate");?>, <?php linkTo("command","horner","horner");?>, <?php linkTo("command","rationalmode","rationalmode");?> 
</div> 
