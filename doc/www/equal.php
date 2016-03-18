<a name="equal"></a> 
<div class="divName"> 
<h2 class="name">Name:</h2> <?php linkTo("command","equal","==");?> 
<span class="smallDescription">equality test operator 
</span> 
</div> 
<div class="divLibraryName"> 
<h2 class="libraryname">Library name:</h2> 
<span class="commandline type">sollya_obj_t sollya_lib_cmp_equal(sollya_obj_t, sollya_obj_t)</span> 
</div> 
<div class="divUsage"> 
<h2 class="category">Usage: </h2> 
<span class="commandline"><span class="arg">expr1</span> <?php linkTo("command","equal","==");?> <span class="arg">expr2</span> : (<span class="type">any type</span>, <span class="type">any type</span>) -&gt; <span class="type">boolean</span></span> 
 
</div> 
<div class="divParameters"> 
<h2 class="category">Parameters: </h2> 
<ul> 
<li><span class="arg">expr1</span> and <span class="arg">expr2</span> represent expressions</li> 
</ul> 
</div> 
<div class="divDescription"> 
<h2 class="category">Description: </h2><ul> 
<li>The test <span class="arg">expr1</span> <?php linkTo("command","equal","==");?> <span class="arg">expr2</span> returns <?php linkTo("command","true","true");?> when <span class="arg">expr1</span> and <span class="arg">expr2</span> are 
syntactically equal and different from <?php linkTo("command","error","error");?> and @NaN@. Conversely if <span class="arg">expr1</span> 
and <span class="arg">expr2</span> are objects that are mathematically different and Sollya manages 
to figure it out, the test returns <?php linkTo("command","false","false");?>. In between these two cases, there 
is the grey zone of expressions that are not syntactically equal but are 
mathematically equal. In such a case, Sollya normally tries to determine if 
the expressions are mathematically equal and if it manages to prove it, it 
returns <?php linkTo("command","true","true");?>, without a warning. In the case when <span class="arg">expr1</span> and <span class="arg">expr2</span> are 
two constant expressions, Sollya will in particular try to evaluate their 
difference: in the case when the difference is 0 or is so small that Sollya 
does not manage to obtain a faithful rounding of the real value, it will 
return <?php linkTo("command","true","true");?> (with a warning if it has not been possible to actually prove 
that the real value is 0). In any other case, when both expressions are not 
syntactically equal and Sollya has not been able to prove that they are 
mathematically equal, it returns <?php linkTo("command","false","false");?>. 
</li><li>The level of simplifications performed by Sollya to determine if 
expressions are mathematically equal depends on the value of <?php linkTo("command","autosimplify","autosimplify");?>. 
If it is <?php linkTo("command","off","off");?>, no formal simplification is performed, hence expression trees 
as simple as x+1 and 1+x will be considered not equal. Conversely, if 
<?php linkTo("command","autosimplify","autosimplify");?> is set to <?php linkTo("command","on","on");?>, polynomial subexpressions that are mathematically 
equal will in general be recognized as being equal. 
</li><li>The user should always keep in mind that a litteral constant written in 
decimal arithmetic (such as 0.1 for instance) is not considered as an exact 
constant by Sollya (unless it is exactly representable in binary without 
requiring too much precision) and is first correctly rounded at precision 
<?php linkTo("command","prec","prec");?>, prior to any other operation. Of course, this leads to a rounding 
warning, but it is important to remember that this is done before the 
expression trees are compared, possibly leading to two expressions comparing 
equal, while they are obviously mathematically different, just because they 
contain different constants that have been rounded to the same value at 
precision <?php linkTo("command","prec","prec");?>. As a general rule, to avoid this behavior, the user should 
represent constants in an exact format such as hexadecimal or represent 
decimal constants as integer fractions (e.g., 0.1 represented by the constant 
expression 1/10). 
</li><li>Notice that @NaN@ and <?php linkTo("command","error","error");?> share the property that they both compare equal 
and different to anything, i.e., if the variable <span class="arg">a</span> contains @NaN@ or <?php linkTo("command","error","error");?> 
and whatever the content of variable <span class="arg">b</span> is, the tests <span class="arg">a</span> <?php linkTo("command","equal","==");?> <span class="arg">b</span> and 
<span class="arg">a</span> <?php linkTo("command","neq","!=");?> <span class="arg">b</span> both return <?php linkTo("command","false","false");?>. The standard way of testing if <span class="arg">a</span> contains 
@NaN@ or <?php linkTo("command","error","error");?> is indeed to check if <span class="arg">a</span> <?php linkTo("command","equal","==");?> <span class="arg">a</span> returns false. In such a 
case, it is however impossible to determine what is the actual value of <span class="arg">a</span> 
amongst both possibilities using only <?php linkTo("command","equal","==");?> or <?php linkTo("command","neq","!=");?>. The standard way to 
discriminate this situation is to use the match ... with ... construct. 
</ul> 
</div> 
<div class="divExamples"> 
<div class="divExample"> 
<h2 class="category">Example 1: </h2> 
&nbsp;&nbsp;&nbsp;&gt; "Hello" == "Hello";<br> 
&nbsp;&nbsp;&nbsp;true<br> 
&nbsp;&nbsp;&nbsp;&gt; "Hello" == "Salut";<br> 
&nbsp;&nbsp;&nbsp;false<br> 
&nbsp;&nbsp;&nbsp;&gt; "Hello" == 5;<br> 
&nbsp;&nbsp;&nbsp;false<br> 
&nbsp;&nbsp;&nbsp;&gt; 5 + x == 5 + x;<br> 
&nbsp;&nbsp;&nbsp;true<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 2: </h2> 
&nbsp;&nbsp;&nbsp;&gt; verbosity = 1!;<br> 
&nbsp;&nbsp;&nbsp;&gt; asin(1) * 2 == pi;<br> 
&nbsp;&nbsp;&nbsp;true<br> 
&nbsp;&nbsp;&nbsp;&gt; cos(3)^2 == 1 - sin(3)^2;<br> 
&nbsp;&nbsp;&nbsp;Warning: the tool is unable to decide an equality test by evaluation even though faithful evaluation of the terms has been possible. The terms will be considered to be equal.<br> 
&nbsp;&nbsp;&nbsp;true<br> 
&nbsp;&nbsp;&nbsp;&gt; exp(5) == log(4);<br> 
&nbsp;&nbsp;&nbsp;false<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 3: </h2> 
&nbsp;&nbsp;&nbsp;&gt; autosimplify=off;<br> 
&nbsp;&nbsp;&nbsp;Automatic pure tree simplification has been deactivated.<br> 
&nbsp;&nbsp;&nbsp;&gt; exp(1+x) == exp(x+1);<br> 
&nbsp;&nbsp;&nbsp;false<br> 
&nbsp;&nbsp;&nbsp;&gt; autosimplify=on;<br> 
&nbsp;&nbsp;&nbsp;Automatic pure tree simplification has been activated.<br> 
&nbsp;&nbsp;&nbsp;&gt; exp(1+x) == exp(x+1);<br> 
&nbsp;&nbsp;&nbsp;false<br> 
&nbsp;&nbsp;&nbsp;&gt; (1/3+x)^2 == x^2 + 1/9 + (5-3)*x/3;<br> 
&nbsp;&nbsp;&nbsp;true<br> 
&nbsp;&nbsp;&nbsp;&gt; log(x)/log(10) == log10(x);<br> 
&nbsp;&nbsp;&nbsp;false<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 4: </h2> 
&nbsp;&nbsp;&nbsp;&gt; prec = 12;<br> 
&nbsp;&nbsp;&nbsp;The precision has been set to 12 bits.<br> 
&nbsp;&nbsp;&nbsp;&gt; verbosity = 1!;<br> 
&nbsp;&nbsp;&nbsp;&gt; 16384.1 == 16385.1;<br> 
&nbsp;&nbsp;&nbsp;Warning: Rounding occurred when converting the constant "16384.1" to floating-point with 12 bits.<br> 
&nbsp;&nbsp;&nbsp;If safe computation is needed, try to increase the precision.<br> 
&nbsp;&nbsp;&nbsp;Warning: Rounding occurred when converting the constant "16385.1" to floating-point with 12 bits.<br> 
&nbsp;&nbsp;&nbsp;If safe computation is needed, try to increase the precision.<br> 
&nbsp;&nbsp;&nbsp;true<br> 
&nbsp;&nbsp;&nbsp;&gt; 16384 == 16384.25;<br> 
&nbsp;&nbsp;&nbsp;false<br> 
&nbsp;&nbsp;&nbsp;&gt; 0.1 == 1/10;<br> 
&nbsp;&nbsp;&nbsp;Warning: Rounding occurred when converting the constant "0.1" to floating-point with 12 bits.<br> 
&nbsp;&nbsp;&nbsp;If safe computation is needed, try to increase the precision.<br> 
&nbsp;&nbsp;&nbsp;false<br> 
&nbsp;&nbsp;&nbsp;&gt; 0.1 == round(1/10, prec, RN);<br> 
&nbsp;&nbsp;&nbsp;Warning: Rounding occurred when converting the constant "0.1" to floating-point with 12 bits.<br> 
&nbsp;&nbsp;&nbsp;If safe computation is needed, try to increase the precision.<br> 
&nbsp;&nbsp;&nbsp;true<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 5: </h2> 
&nbsp;&nbsp;&nbsp;&gt; error == error;<br> 
&nbsp;&nbsp;&nbsp;false<br> 
&nbsp;&nbsp;&nbsp;&gt; error != error;<br> 
&nbsp;&nbsp;&nbsp;false<br> 
&nbsp;&nbsp;&nbsp;&gt; @NaN@ == @NaN@;<br> 
&nbsp;&nbsp;&nbsp;false<br> 
&nbsp;&nbsp;&nbsp;&gt; @NaN@ != @NaN@;<br> 
&nbsp;&nbsp;&nbsp;false<br> 
&nbsp;&nbsp;&nbsp;&gt; error == @NaN@;<br> 
&nbsp;&nbsp;&nbsp;false<br> 
&nbsp;&nbsp;&nbsp;&gt; error != @NaN@;<br> 
&nbsp;&nbsp;&nbsp;false<br> 
&nbsp;&nbsp;&nbsp;&gt; a = error;<br> 
&nbsp;&nbsp;&nbsp;&gt; match a with<br> 
&nbsp;&nbsp;&nbsp;&nbsp;  @NaN@ : ("a contains @NaN@")<br> 
&nbsp;&nbsp;&nbsp;&nbsp;  default:("a contains something else");<br> 
&nbsp;&nbsp;&nbsp;error<br> 
&nbsp;&nbsp;&nbsp;&gt; a = @NaN@;<br> 
&nbsp;&nbsp;&nbsp;&gt; match a with<br> 
&nbsp;&nbsp;&nbsp;&nbsp;  @NaN@ : ("a contains @NaN@")<br> 
&nbsp;&nbsp;&nbsp;&nbsp;  default:("a contains something else");<br> 
&nbsp;&nbsp;&nbsp;a contains @NaN@<br> 
</div> 
</div> 
<div class="divSeeAlso"> 
<span class="category">See also: </span><?php linkTo("command","neq","!=");?>, <?php linkTo("command","gt","&gt;");?>, <?php linkTo("command","ge","&gt;=");?>, <?php linkTo("command","le","&lt;=");?>, <?php linkTo("command","lt","&lt;");?>, <?php linkTo("command","in","in");?>, <?php linkTo("command","not","!");?>, <?php linkTo("command","and","&&");?>, <?php linkTo("command","or","||");?>, <?php linkTo("command","error","error");?>, <?php linkTo("command","prec","prec");?>, <?php linkTo("command","autosimplify","autosimplify");?> 
</div> 
