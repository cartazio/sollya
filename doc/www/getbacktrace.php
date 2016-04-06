<a name="getbacktrace"></a> 
<div class="divName"> 
<h2 class="name">Name:</h2> <?php linkTo("command","getbacktrace","getbacktrace");?> 
<span class="smallDescription">returns the list of Sollya procedures currently run 
</span> 
</div> 
<div class="divLibraryName"> 
<h2 class="libraryname">Library name:</h2> 
<span class="commandline type">sollya_obj_t sollya_lib_getbacktrace();</span> 
</div> 
<div class="divUsage"> 
<h2 class="category">Usage: </h2> 
<span class="commandline"><?php linkTo("command","getbacktrace","getbacktrace");?>() : <span class="type">void</span> -&gt; <span class="type">list</span></span> 
 
</div> 
<div class="divDescription"> 
<h2 class="category">Description: </h2><ul> 
<li>The <?php linkTo("command","getbacktrace","getbacktrace");?> command allows the stack of Sollya procedures that are 
currently run to be inspected. When called, <?php linkTo("command","getbacktrace","getbacktrace");?> returns an 
ordered list of structures, each of which contains an element 
passed_args and an element called_proc. The element called_proc 
contains the Sollya object representing the procedure being run. The 
element passed_args contains an ordered list of all effective 
arguments passed to the procedure when it was called. The procedure called 
last (i.e., on top of the stack) comes first in the list returned 
by <?php linkTo("command","getbacktrace","getbacktrace");?>. When any of the procedure called takes no arguments, the 
passed_arg element of the corresponding structure evaluates to an empty 
list. 
</li><li>When called from outside any procedure (at toplevel), <?php linkTo("command","getbacktrace","getbacktrace");?> returns 
an empty list. 
</li><li>When called for a stack containing a call to a variadic procedure that was 
called with an infinite number of effective arguments, the corresponding 
passed_args element evaluates to an end-elliptic list. 
</ul> 
</div> 
<div class="divExamples"> 
<div class="divExample"> 
<h2 class="category">Example 1: </h2> 
&nbsp;&nbsp;&nbsp;&gt; procedure testA() {<br> 
&nbsp;&nbsp;&nbsp;&nbsp; 	&nbsp;&nbsp;"Current backtrace:";<br> 
&nbsp;&nbsp;&nbsp;&nbsp; 	&nbsp;&nbsp;getbacktrace();<br> 
&nbsp;&nbsp;&nbsp;&nbsp; };<br> 
&nbsp;&nbsp;&nbsp;&gt; procedure testB(X) {<br> 
&nbsp;&nbsp;&nbsp;&nbsp; 	&nbsp;&nbsp;"X = ", X;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; 	&nbsp;&nbsp;testA();<br> 
&nbsp;&nbsp;&nbsp;&nbsp; };<br> 
&nbsp;&nbsp;&nbsp;&gt; procedure testC(X, Y) {<br> 
&nbsp;&nbsp;&nbsp;&nbsp; 	&nbsp;&nbsp;"X = ", X, ", Y = ", Y;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; 	&nbsp;&nbsp;testB(Y);<br> 
&nbsp;&nbsp;&nbsp;&nbsp; };<br> 
&nbsp;&nbsp;&nbsp;&gt; testC(17, 42);<br> 
&nbsp;&nbsp;&nbsp;X = 17, Y = 42<br> 
&nbsp;&nbsp;&nbsp;X = 42<br> 
&nbsp;&nbsp;&nbsp;Current backtrace:<br> 
&nbsp;&nbsp;&nbsp;[|{ .passed_args = [| |], .called_proc = proc()<br> 
&nbsp;&nbsp;&nbsp;{<br> 
&nbsp;&nbsp;&nbsp;"Current backtrace:";<br> 
&nbsp;&nbsp;&nbsp;getbacktrace();<br> 
&nbsp;&nbsp;&nbsp;return void;<br> 
&nbsp;&nbsp;&nbsp;} }, { .passed_args = [|42|], .called_proc = proc(X)<br> 
&nbsp;&nbsp;&nbsp;{<br> 
&nbsp;&nbsp;&nbsp;"X = ", X;<br> 
&nbsp;&nbsp;&nbsp;testA();<br> 
&nbsp;&nbsp;&nbsp;return void;<br> 
&nbsp;&nbsp;&nbsp;} }, { .passed_args = [|17, 42|], .called_proc = proc(X, Y)<br> 
&nbsp;&nbsp;&nbsp;{<br> 
&nbsp;&nbsp;&nbsp;"X = ", X, ", Y = ", Y;<br> 
&nbsp;&nbsp;&nbsp;testB(Y);<br> 
&nbsp;&nbsp;&nbsp;return void;<br> 
&nbsp;&nbsp;&nbsp;} }|]<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 2: </h2> 
&nbsp;&nbsp;&nbsp;&gt; getbacktrace();<br> 
&nbsp;&nbsp;&nbsp;[| |]<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 3: </h2> 
&nbsp;&nbsp;&nbsp;&gt; procedure printnumargs(X) {<br> 
&nbsp;&nbsp;&nbsp;&nbsp; 	&nbsp;&nbsp;var L, t;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; 	&nbsp;&nbsp;"number of arguments: ", X;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; 	&nbsp;&nbsp;L = getbacktrace();<br> 
&nbsp;&nbsp;&nbsp;&nbsp; 	&nbsp;&nbsp;"Backtrace:";<br> 
&nbsp;&nbsp;&nbsp;&nbsp; 	&nbsp;&nbsp;for t in L do {<br> 
&nbsp;&nbsp;&nbsp;&nbsp; 	&nbsp;&nbsp;&nbsp;&nbsp;"&nbsp;&nbsp;" @ objectname(t.called_proc) @ ", ", t.passed_args;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; 	&nbsp;&nbsp;};<br> 
&nbsp;&nbsp;&nbsp;&nbsp; };<br> 
&nbsp;&nbsp;&nbsp;&gt; procedure numargs(l = ...) {<br> 
&nbsp;&nbsp;&nbsp;&nbsp; 	&nbsp;&nbsp;"l[17] = ", l[17];<br> 
&nbsp;&nbsp;&nbsp;&nbsp; 	&nbsp;&nbsp;printnumargs(length(l));<br> 
&nbsp;&nbsp;&nbsp;&nbsp; };<br> 
&nbsp;&nbsp;&nbsp;&gt; procedure test() {<br> 
&nbsp;&nbsp;&nbsp;&nbsp; 	&nbsp;&nbsp;numargs @ [|25, 26, 27 ...|];<br> 
&nbsp;&nbsp;&nbsp;&nbsp; };<br> 
&nbsp;&nbsp;&nbsp;&gt; test();<br> 
&nbsp;&nbsp;&nbsp;l[17] = 42<br> 
&nbsp;&nbsp;&nbsp;number of arguments: infty<br> 
&nbsp;&nbsp;&nbsp;Backtrace:<br> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;printnumargs, [|infty|]<br> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;numargs, [|25, 26, 27...|]<br> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;test, [| |]<br> 
</div> 
</div> 
<div class="divSeeAlso"> 
<span class="category">See also: </span><?php linkTo("command","proc","proc");?>, <?php linkTo("command","procedure","procedure");?>, <?php linkTo("command","objectname","objectname");?>, <?php linkTo("command","bind","bind");?>, <?php linkTo("command","concat","@");?> 
</div> 
