<div class="divExample">
&nbsp;&nbsp;&nbsp;&gt; f = exp(x);<br>
&nbsp;&nbsp;&nbsp;&gt; f;<br>
&nbsp;&nbsp;&nbsp;exp(x)<br>
&nbsp;&nbsp;&nbsp;&gt; a = "Hello world";<br>
&nbsp;&nbsp;&nbsp;&gt; a;<br>
&nbsp;&nbsp;&nbsp;Hello world<br>
&nbsp;&nbsp;&nbsp;&gt; b = 5;<br>
&nbsp;&nbsp;&nbsp;&gt; f(b);<br>
&nbsp;&nbsp;&nbsp;Warning: rounding has happened. The value displayed is a faithful rounding to 165 bits of the true result.<br>
&nbsp;&nbsp;&nbsp;148.41315910257660342111558004055227962348766759388<br>
&nbsp;&nbsp;&nbsp;&gt; {var b; b = 4; f(b); };<br>
&nbsp;&nbsp;&nbsp;Warning: rounding has happened. The value displayed is a faithful rounding to 165 bits of the true result.<br>
&nbsp;&nbsp;&nbsp;54.598150033144239078110261202860878402790737038614<br>
&nbsp;&nbsp;&nbsp;&gt; {var x; x = 3; };<br>
&nbsp;&nbsp;&nbsp;Warning: the identifier "x" is already bound to the current free variable.<br>
&nbsp;&nbsp;&nbsp;It cannot be declared as a local variable. The declaration of "x" will have no effect.<br>
&nbsp;&nbsp;&nbsp;Warning: the identifier "x" is already bound to the free variable, to a library function, library constant or to an external procedure.<br>
&nbsp;&nbsp;&nbsp;The command will have no effect.<br>
&nbsp;&nbsp;&nbsp;Warning: the last assignment will have no effect.<br>
&nbsp;&nbsp;&nbsp;&gt; {var a, b; a=5; b=3; {var a; var b; b = true; a = 1; a; b;}; a; b; };<br>
&nbsp;&nbsp;&nbsp;1<br>
&nbsp;&nbsp;&nbsp;true<br>
&nbsp;&nbsp;&nbsp;5<br>
&nbsp;&nbsp;&nbsp;3<br>
&nbsp;&nbsp;&nbsp;&gt; a;<br>
&nbsp;&nbsp;&nbsp;Hello world<br>
</div>
