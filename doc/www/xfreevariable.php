<a name="xfreevariable"></a> 
<div class="divName"> 
<h2 class="name">Name:</h2> <?php linkTo("command","xfreevariable","_x_");?> 
<span class="smallDescription">universal name for the mathematical free variable. 
</span> 
</div> 
<div class="divLibraryName"> 
<h2 class="libraryname">Library names:</h2> 
<span class="commandline type">sollya_obj_t sollya_lib_free_variable()</span> 
<span class="commandline type">sollya_obj_t sollya_lib_build_function_free_variable()</span> 
<span class="commandline type">#define SOLLYA_X_ (sollya_lib_build_function_free_variable())</span> 
</div> 
<div class="divDescription"> 
<h2 class="category">Description: </h2><ul> 
<li><?php linkTo("command","xfreevariable","_x_");?> is an identifier that always denotes the mathematical free variable. 
It cannot be assigned. 
</li><li>Sollya manipulates mathematical functions of a single variable. The first 
time that a variable name is used without having been assigned before, this 
variable name is automatically considered by Sollya as the name of the 
free variable. Subsequently, any other unassigned variable name will be 
considered as the free variable with a warning making this conversion 
explicit. This is convenient for an every-day use of the interactive tool, 
but it has the drawback that the free variable name can change from a 
session to another. There are contexts (e.g., within a procedure, or for 
doing pattern matching) when one might want to refer to the free variable 
regardless of its name in the current session. For this purpose <?php linkTo("command","xfreevariable","_x_");?> is 
a universal identifier, always available and always denoting the free 
variable, whatever its name is in the current context. 
</ul> 
</div> 
<div class="divExamples"> 
<div class="divExample"> 
<h2 class="category">Example 1: </h2> 
&nbsp;&nbsp;&nbsp;&gt; verbosity=1!;<br> 
&nbsp;&nbsp;&nbsp;&gt; sin(a);<br> 
&nbsp;&nbsp;&nbsp;sin(a)<br> 
&nbsp;&nbsp;&nbsp;&gt; b;<br> 
&nbsp;&nbsp;&nbsp;Warning: the identifier "b" is neither assigned to, nor bound to a library function nor external procedure, nor equal to the current free variable.<br> 
&nbsp;&nbsp;&nbsp;Will interpret "b" as "a".<br> 
&nbsp;&nbsp;&nbsp;a<br> 
&nbsp;&nbsp;&nbsp;&gt; _x_;<br> 
&nbsp;&nbsp;&nbsp;a<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 2: </h2> 
&nbsp;&nbsp;&nbsp;&gt; verbosity=1!;<br> 
&nbsp;&nbsp;&nbsp;&gt; sin(y);<br> 
&nbsp;&nbsp;&nbsp;sin(y)<br> 
&nbsp;&nbsp;&nbsp;&gt; f = proc(a) {<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;return sin(a + _x_);<br> 
&nbsp;&nbsp;&nbsp;&nbsp; };<br> 
&nbsp;&nbsp;&nbsp;&gt; rename(y,z);<br> 
&nbsp;&nbsp;&nbsp;Information: the free variable has been renamed from "y" to "z".<br> 
&nbsp;&nbsp;&nbsp;&gt; f(1);<br> 
&nbsp;&nbsp;&nbsp;sin(1 + z)<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 3: </h2> 
&nbsp;&nbsp;&nbsp;&gt; f = sin(y);<br> 
&nbsp;&nbsp;&nbsp;&gt; match f with<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;sin(a) : { print("sin of a with a =", a);<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;match a with<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;_x_ : { print("a turns out to be the free variable"); }<br> 
&nbsp;&nbsp;&nbsp;&nbsp; 	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;default : { print("a is some expression"); };<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;}<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;_x_ : { print("Free variable") ; }<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;default: { print("Something else"); };<br> 
&nbsp;&nbsp;&nbsp;sin of a with a = y<br> 
&nbsp;&nbsp;&nbsp;a turns out to be the free variable<br> 
</div> 
</div> 
<div class="divSeeAlso"> 
<span class="category">See also: </span><?php linkTo("command","rename","rename");?>, <?php linkTo("command","isbound","isbound");?>, <?php linkTo("command","proc","proc");?> 
</div> 
