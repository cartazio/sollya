<a name="objectname"></a> 
<div class="divName"> 
<h2 class="name">Name:</h2> <?php linkTo("command","objectname","objectname");?> 
<span class="smallDescription">returns, given a Sollya object, a string that can be reparsed to the object 
</span> 
</div> 
<div class="divLibraryName"> 
<h2 class="libraryname">Library name:</h2> 
<span class="commandline type">sollya_obj_t sollya_lib_objectname(sollya_obj_t);</span> 
</div> 
<div class="divUsage"> 
<h2 class="category">Usage: </h2> 
<span class="commandline"><?php linkTo("command","objectname","objectname");?>(<span class="arg">obj</span>) : <span class="type">any type</span> -&gt; <span class="type">string</span></span> 
 
</div> 
<div class="divDescription"> 
<h2 class="category">Description: </h2><ul> 
<li><?php linkTo("command","objectname","objectname");?>(<span class="arg">obj</span>) queries the Sollya symbol table in order to recover the 
name of an identifier the object <span class="arg">obj</span> is assigned to. If it succeeds, it 
returns a string containing the recovered identifier. In contrast, if it 
does not succeed, it returns a string simply containing a textual 
representation of <span class="arg">obj</span>. 
</li><li>The only condition for an identifier to be eligible to be returned by 
<?php linkTo("command","objectname","objectname");?>(<span class="arg">obj</span>) is to be accessible in the scope <?php linkTo("command","objectname","objectname");?> is executed in, 
i.e., not to be shadowed by an identifier of the same name which does not 
hold the object <span class="arg">obj</span>. 
</li><li>In any case, if the string returned by <?php linkTo("command","objectname","objectname");?> is given to the <?php linkTo("command","parse","parse");?> 
command in the same scope, the original object <span class="arg">obj</span> is recovered. 
</li><li><?php linkTo("command","objectname","objectname");?> is particularly useful in combination with <?php linkTo("command","getbacktrace","getbacktrace");?>, when 
the Sollya procedure stack is to be displayed in a fashion, where 
procedures are identified by their name and not their procedural content. 
</li><li><?php linkTo("command","objectname","objectname");?> may also be used to get a string representation of the free 
mathematical variable. 
</li><li>If an object is simply to be cast into a string, without trying to 
retrieve an identifier for it, <?php linkTo("command","objectname","objectname");?> is not appropriate. In this case, 
it suffices to concatenate it to an empty string with the <?php linkTo("command","concat","@");?> operator. 
</ul> 
</div> 
<div class="divExamples"> 
<div class="divExample"> 
<h2 class="category">Example 1: </h2> 
&nbsp;&nbsp;&nbsp;&gt; s = "Hello";<br> 
&nbsp;&nbsp;&nbsp;&gt; objectname("Hello");<br> 
&nbsp;&nbsp;&nbsp;s<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 2: </h2> 
&nbsp;&nbsp;&nbsp;&gt; f = exp(x);<br> 
&nbsp;&nbsp;&nbsp;&gt; g = sin(x);<br> 
&nbsp;&nbsp;&nbsp;&gt; [| objectname(exp(x)), objectname(sin(x)), objectname(cos(x)) |];<br> 
&nbsp;&nbsp;&nbsp;[|"f", "g", "cos(x)"|]<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 3: </h2> 
&nbsp;&nbsp;&nbsp;&gt; o = { .f = exp(x), .I = [-1;1] };<br> 
&nbsp;&nbsp;&nbsp;&gt; s1 = o@""; s1;<br> 
&nbsp;&nbsp;&nbsp;{ .f = exp(x), .I = [-1;1] }<br> 
&nbsp;&nbsp;&nbsp;&gt; s2 = objectname({ .I = [-1;1], .f = exp(x)}); s2;<br> 
&nbsp;&nbsp;&nbsp;o<br> 
&nbsp;&nbsp;&nbsp;&gt; parse(s1) == parse(s2);<br> 
&nbsp;&nbsp;&nbsp;true<br> 
&nbsp;&nbsp;&nbsp;&gt; write("s2 = \"", s2, "\" parses to ", parse(s2), "\n");<br> 
&nbsp;&nbsp;&nbsp;s2 = "o" parses to { .f = exp(x), .I = [-1;1] }<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 4: </h2> 
&nbsp;&nbsp;&nbsp;&gt; n = 1664;<br> 
&nbsp;&nbsp;&nbsp;&gt; objectname(n);<br> 
&nbsp;&nbsp;&nbsp;n<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 5: </h2> 
&nbsp;&nbsp;&nbsp;&gt; f = exp(x);<br> 
&nbsp;&nbsp;&nbsp;&gt; g = sin(x);<br> 
&nbsp;&nbsp;&nbsp;&gt; procedure test() {<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;var f;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;var h;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;f = tan(x);<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;h = cos(x);<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;[| objectname(exp(x)), objectname(sin(x)), objectname(cos(x)), objectname(tan(x)) |];<br> 
&nbsp;&nbsp;&nbsp;&nbsp; };<br> 
&nbsp;&nbsp;&nbsp;&gt; test();<br> 
&nbsp;&nbsp;&nbsp;[|"exp(x)", "g", "h", "f"|]<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 6: </h2> 
&nbsp;&nbsp;&nbsp;&gt; procedure apply_proc(p, a, b) {<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;return p(a, b);<br> 
&nbsp;&nbsp;&nbsp;&nbsp; };<br> 
&nbsp;&nbsp;&nbsp;&gt; procedure show_trace_and_add(n, m) {<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;var i, bt;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;bt = getbacktrace();<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;write("Procedure stack:\n");<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;for i from 0 to length(bt) - 1 do {<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;	write("&nbsp;&nbsp;&nbsp;Procedure ", objectname((bt[i]).called_proc), " called with ", length((bt[i]).passed_args), " arguments\n");<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;};<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;write("\n");<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;return n + m;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; };<br> 
&nbsp;&nbsp;&nbsp;&gt; procedure show_and_succ(u) {<br> 
&nbsp;&nbsp;&nbsp;&nbsp; 	&nbsp;&nbsp;return apply_proc(show_trace_and_add, u, 1);<br> 
&nbsp;&nbsp;&nbsp;&nbsp; };<br> 
&nbsp;&nbsp;&nbsp;&gt; show_and_succ(16);<br> 
&nbsp;&nbsp;&nbsp;Procedure stack:<br> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Procedure show_trace_and_add called with 2 arguments<br> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Procedure apply_proc called with 3 arguments<br> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Procedure show_and_succ called with 1 arguments<br> 
&nbsp;&nbsp;&nbsp;<br> 
&nbsp;&nbsp;&nbsp;17<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 7: </h2> 
&nbsp;&nbsp;&nbsp;&gt; f = exp(three_decker_sauerkraut_and_toadstool_sandwich_with_arsenic_sauce);<br> 
&nbsp;&nbsp;&nbsp;&gt; g = sin(_x_);<br> 
&nbsp;&nbsp;&nbsp;&gt; h = f(g);<br> 
&nbsp;&nbsp;&nbsp;&gt; h;<br> 
&nbsp;&nbsp;&nbsp;exp(sin(three_decker_sauerkraut_and_toadstool_sandwich_with_arsenic_sauce))<br> 
&nbsp;&nbsp;&nbsp;&gt; objectname(_x_);<br> 
&nbsp;&nbsp;&nbsp;three_decker_sauerkraut_and_toadstool_sandwich_with_arsenic_sauce<br> 
</div> 
</div> 
<div class="divSeeAlso"> 
<span class="category">See also: </span><?php linkTo("command","parse","parse");?>, <?php linkTo("command","var","var");?>, <?php linkTo("command","getbacktrace","getbacktrace");?>, <?php linkTo("command","proc","proc");?>, <?php linkTo("command","procedure","procedure");?>, <?php linkTo("command","concat","@");?> 
</div> 
