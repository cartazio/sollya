<div class="divExample">
&nbsp;&nbsp;&nbsp;&gt; f = exp(x) + 5;<br>
&nbsp;&nbsp;&nbsp;&gt; f(-infty);<br>
&nbsp;&nbsp;&nbsp;5<br>
&nbsp;&nbsp;&nbsp;&gt; evaluate(f,[-infty;infty]);<br>
&nbsp;&nbsp;&nbsp;[5;infty]<br>
&nbsp;&nbsp;&nbsp;&gt; f(infty);<br>
&nbsp;&nbsp;&nbsp;infty<br>
&nbsp;&nbsp;&nbsp;&gt; [-infty;5] * [3;4];<br>
&nbsp;&nbsp;&nbsp;[-infty;20]<br>
&nbsp;&nbsp;&nbsp;&gt; -infty < 5;<br>
&nbsp;&nbsp;&nbsp;true<br>
&nbsp;&nbsp;&nbsp;&gt; log(0);<br>
&nbsp;&nbsp;&nbsp;-infty<br>
&nbsp;&nbsp;&nbsp;&gt; [log(0);17];<br>
&nbsp;&nbsp;&nbsp;Warning: inclusion property is satisfied but the diameter may be greater than the least possible.<br>
&nbsp;&nbsp;&nbsp;Warning: at least one of the given expressions is not a constant but requires evaluation.<br>
&nbsp;&nbsp;&nbsp;Evaluation is guaranteed to ensure the inclusion property. The approximate result is at least 165 bit accurate.<br>
&nbsp;&nbsp;&nbsp;[-infty;17]<br>
&nbsp;&nbsp;&nbsp;&gt; <br>
</div>
