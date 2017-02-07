<a name="eucldiv"></a> 
<div class="divName"> 
<h2 class="name">Name:</h2> <?php linkTo("command","eucldiv","div");?> 
<span class="smallDescription">Computes the euclidian division of polynomials or numbers and returns the quotient 
</span> 
</div> 
<div class="divLibraryName"> 
<h2 class="libraryname">Library name:</h2> 
<span class="commandline type">sollya_obj_t sollya_lib_euclidian_div(sollya_obj_t, sollya_obj_t)</span> 
</div> 
<div class="divUsage"> 
<h2 class="category">Usage: </h2> 
<span class="commandline"><?php linkTo("command","eucldiv","div");?>(<span class="arg">a</span>, <span class="arg">b</span>) : (<span class="type">function</span>, <span class="type">function</span>) -&gt; <span class="type">function</span></span> 
 
</div> 
<div class="divParameters"> 
<h2 class="category">Parameters: </h2> 
<ul> 
<li><span class="arg">a</span> is a polynomial.</li> 
<li><span class="arg">b</span> is a polynomial.</li> 
</ul> 
</div> 
<div class="divDescription"> 
<h2 class="category">Description: </h2><ul> 
<li>When both <span class="arg">a</span> and <span class="arg">b</span> are constants, <?php linkTo("command","eucldiv","div");?>(<span class="arg">a</span>,<span class="arg">b</span>) computes the 
largest integer less than or equal to <span class="arg">a</span> divided by <span class="arg">b</span>. In other 
words, it returns the quotient of the Euclidian division of <span class="arg">a</span> by 
<span class="arg">b</span>. 
</li><li>When at least one of <span class="arg">a</span> or <span class="arg">b</span> is a polynomial of degree at least 
1, <?php linkTo("command","eucldiv","div");?>(<span class="arg">a</span>,<span class="arg">b</span>) computes two polynomials <span class="arg">q</span> and <span class="arg">r</span> such 
that <span class="arg">a</span> is equal to the product of <span class="arg">q</span> and <span class="arg">b</span> plus <span class="arg">r</span>. The 
polynomial <span class="arg">r</span> is of least degree possible. The <?php linkTo("command","eucldiv","div");?> command 
returns <span class="arg">q</span>. In order to recover <span class="arg">r</span>, use the <?php linkTo("command","euclmod","mod");?> command. 
</li><li>When at least one of <span class="arg">a</span> or <span class="arg">b</span> is a function that is no polynomial, 
<?php linkTo("command","eucldiv","div");?>(<span class="arg">a</span>,<span class="arg">b</span>) returns 0. 
</ul> 
</div> 
<div class="divExamples"> 
<div class="divExample"> 
<h2 class="category">Example 1: </h2> 
&nbsp;&nbsp;&nbsp;&gt; div(1001, 231);<br> 
&nbsp;&nbsp;&nbsp;4<br> 
&nbsp;&nbsp;&nbsp;&gt; div(13, 17);<br> 
&nbsp;&nbsp;&nbsp;0<br> 
&nbsp;&nbsp;&nbsp;&gt; div(-14, 15);<br> 
&nbsp;&nbsp;&nbsp;-1<br> 
&nbsp;&nbsp;&nbsp;&gt; div(-213, -5);<br> 
&nbsp;&nbsp;&nbsp;42<br> 
&nbsp;&nbsp;&nbsp;&gt; div(23/13, 11/17);<br> 
&nbsp;&nbsp;&nbsp;2<br> 
&nbsp;&nbsp;&nbsp;&gt; div(exp(13),-sin(17));<br> 
&nbsp;&nbsp;&nbsp;460177<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 2: </h2> 
&nbsp;&nbsp;&nbsp;&gt; div(24 + 68 * x + 74 * x^2 + 39 * x^3 + 10 * x^4 + x^5, 4 + 4 * x + x^2);<br> 
&nbsp;&nbsp;&nbsp;6 + x * (11 + x * (6 + x))<br> 
&nbsp;&nbsp;&nbsp;&gt; div(24 + 68 * x + 74 * x^2 + 39 * x^3 + 10 * x^4 + x^5, 2 * x^3);<br> 
&nbsp;&nbsp;&nbsp;19.5 + x * (5 + x * 0.5)<br> 
&nbsp;&nbsp;&nbsp;&gt; div(x^2, x^3);<br> 
&nbsp;&nbsp;&nbsp;0<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 3: </h2> 
&nbsp;&nbsp;&nbsp;&gt; div(exp(x), x^2);<br> 
&nbsp;&nbsp;&nbsp;0<br> 
</div> 
</div> 
<div class="divSeeAlso"> 
<span class="category">See also: </span><?php linkTo("command","gcd","gcd");?>, <?php linkTo("command","euclmod","mod");?>, <?php linkTo("command","numberroots","numberroots");?> 
</div> 
