<div class="divExample">
&nbsp;&nbsp;&nbsp;&gt; evaluate(exp(x),[-infty;0]);<br>
&nbsp;&nbsp;&nbsp;[0;1]<br>
&nbsp;&nbsp;&nbsp;&gt; dirtyinfnorm(exp(x),[-infty;0]);<br>
&nbsp;&nbsp;&nbsp;Warning: a bound of the interval is infinite or NaN.<br>
&nbsp;&nbsp;&nbsp;This command cannot handle such intervals.<br>
&nbsp;&nbsp;&nbsp;NaN<br>
&nbsp;&nbsp;&nbsp;&gt; <br>
&nbsp;&nbsp;&nbsp;&gt; f = log(x);<br>
&nbsp;&nbsp;&nbsp;&gt; [f(0); f(1)];<br>
&nbsp;&nbsp;&nbsp;Warning: inclusion property is satisfied but the diameter may be greater than the least possible.<br>
&nbsp;&nbsp;&nbsp;Warning: at least one of the given expressions is not a constant but requires evaluation.<br>
&nbsp;&nbsp;&nbsp;Evaluation is guaranteed to ensure the inclusion property. The approximate result is at least 165 bit accurate.<br>
&nbsp;&nbsp;&nbsp;[-infty;0]<br>
&nbsp;&nbsp;&nbsp;&gt; <br>
</div>
