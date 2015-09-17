<div class="divExample">
&nbsp;&nbsp;&nbsp;&gt; prec = 300!;<br>
&nbsp;&nbsp;&nbsp;&gt; a = 4097.1;<br>
&nbsp;&nbsp;&nbsp;Warning: Rounding occurred when converting the constant "4097.1" to floating-point with 300 bits.<br>
&nbsp;&nbsp;&nbsp;If safe computation is needed, try to increase the precision.<br>
&nbsp;&nbsp;&nbsp;&gt; prec = 12!;<br>
&nbsp;&nbsp;&nbsp;&gt; d = [4097.1; a];<br>
&nbsp;&nbsp;&nbsp;Warning: Rounding occurred when converting the constant "4097.1" to floating-point with 12 bits.<br>
&nbsp;&nbsp;&nbsp;The inclusion property is nevertheless satisfied.<br>
&nbsp;&nbsp;&nbsp;&gt; prec = 300!;<br>
&nbsp;&nbsp;&nbsp;&gt; d;<br>
&nbsp;&nbsp;&nbsp;[4096;4097.1]<br>
&nbsp;&nbsp;&nbsp;&gt; prec = 30!;<br>
&nbsp;&nbsp;&nbsp;&gt; [-pi;pi];<br>
&nbsp;&nbsp;&nbsp;Warning: at least one of the given expressions is not a constant but requires evaluation.<br>
&nbsp;&nbsp;&nbsp;Evaluation is guaranteed to ensure the inclusion property. The approximate result is at least 30 bit accurate.<br>
&nbsp;&nbsp;&nbsp;[-3.141592655;3.141592655]<br>
</div>
