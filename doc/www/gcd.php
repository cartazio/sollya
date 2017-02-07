<a name="gcd"></a> 
<div class="divName"> 
<h2 class="name">Name:</h2> <?php linkTo("command","gcd","gcd");?> 
<span class="smallDescription">Computes the greatest common divisor of polynomials or numbers. 
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
<li>When both <span class="arg">p</span> and <span class="arg">q</span> are integers, <?php linkTo("command","gcd","gcd");?>(<span class="arg">p</span>,<span class="arg">q</span>) computes the 
greatest common divisor of these two integers, i.e. the greatest 
integer dividing both <span class="arg">p</span> and <span class="arg">q</span>. 
</li><li>When both <span class="arg">p</span> and <span class="arg">q</span> are rational numbers, say a/b and 
c/d, <?php linkTo("command","gcd","gcd");?>(<span class="arg">p</span>,<span class="arg">q</span>) computes the greatest common divisor 
of a * d and b * c, divided by the 
product of the denominators, b * d. 
</li><li>When both <span class="arg">p</span> and <span class="arg">q</span> are constants but at least one of them is 
no rational number, <?php linkTo("command","gcd","gcd");?>(<span class="arg">p</span>,<span class="arg">q</span>) returns 1. 
</li><li>When at least one of <span class="arg">p</span> or <span class="arg">q</span> is a polynomial of degree at least 
1, <?php linkTo("command","gcd","gcd");?>(<span class="arg">p</span>,<span class="arg">q</span>) returns the polynomial of greatest degree 
dividing both <span class="arg">p</span> and <span class="arg">q</span>. The coefficient of the monomial of greatest 
degree of the returned polynomial is set to the greatest common 
divisor of the coefficients of the monomials of greatest degree of <span class="arg">p</span> 
and <span class="arg">q</span>. 
</li><li>When at least one of <span class="arg">p</span> or <span class="arg">q</span> is a function that is no polynomial, 
<?php linkTo("command","gcd","gcd");?>(<span class="arg">p</span>,<span class="arg">q</span>) returns 1. 
</ul> 
</div> 
<div class="divExamples"> 
<div class="divExample"> 
<h2 class="category">Example 1: </h2> 
&nbsp;&nbsp;&nbsp;&gt; gcd(1001, 231);<br> 
&nbsp;&nbsp;&nbsp;77<br> 
&nbsp;&nbsp;&nbsp;&gt; gcd(13, 17);<br> 
&nbsp;&nbsp;&nbsp;1<br> 
&nbsp;&nbsp;&nbsp;&gt; gcd(-210, 462);<br> 
&nbsp;&nbsp;&nbsp;42<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 2: </h2> 
&nbsp;&nbsp;&nbsp;&gt; rationalmode = on!;<br> 
&nbsp;&nbsp;&nbsp;&gt; gcd(6/7, 33/13);<br> 
&nbsp;&nbsp;&nbsp;3 / 91<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 3: </h2> 
&nbsp;&nbsp;&nbsp;&gt; gcd(exp(13),sin(17));<br> 
&nbsp;&nbsp;&nbsp;1<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 4: </h2> 
&nbsp;&nbsp;&nbsp;&gt; gcd(24 + 68 * x + 74 * x^2 + 39 * x^3 + 10 * x^4 + x^5, 480 + 776 * x + 476 * x^2 + 138 * x^3 + 19 * x^4 + x^5);<br> 
&nbsp;&nbsp;&nbsp;4 + x * (4 + x)<br> 
&nbsp;&nbsp;&nbsp;&gt; gcd(1001 * x^2, 231 * x);<br> 
&nbsp;&nbsp;&nbsp;x * 77<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 5: </h2> 
&nbsp;&nbsp;&nbsp;&gt; gcd(exp(x), x^2);<br> 
&nbsp;&nbsp;&nbsp;1<br> 
</div> 
</div> 
<div class="divSeeAlso"> 
<span class="category">See also: </span><?php linkTo("command","eucldiv","div");?>, <?php linkTo("command","euclmod","mod");?>, <?php linkTo("command","numberroots","numberroots");?> 
</div> 
