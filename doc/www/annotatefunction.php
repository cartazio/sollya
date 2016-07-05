<a name="annotatefunction"></a> 
<div class="divName"> 
<h2 class="name">Name:</h2> <?php linkTo("command","annotatefunction","annotatefunction");?> 
<span class="smallDescription">Annotates a Sollya function object with an approximation that is faster to evaluate 
</span> 
</div> 
<div class="divLibraryName"> 
<h2 class="libraryname">Library names:</h2> 
<span class="commandline type">sollya_obj_t sollya_lib_annotatefunction(sollya_obj_t, sollya_obj_t,</span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sollya_obj_t, sollya_obj_t, ...);</span> 
<span class="commandline type">sollya_obj_t sollya_lib_v_annotatefunction(sollya_obj_t, sollya_obj_t,</span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sollya_obj_t, sollya_obj_t,</span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;va_list);</span> 
</div> 
<div class="divUsage"> 
<h2 class="category">Usage: </h2> 
<span class="commandline"><?php linkTo("command","annotatefunction","annotatefunction");?>(<span class="arg">f</span>, <span class="arg">g</span>, <span class="arg">I</span>, <span class="arg">d</span>) : (<span class="type">function</span>, <span class="type">function</span>, <span class="type">range</span>, <span class="type">range</span>) -&gt; <span class="type">function</span></span> 
<span class="commandline"><?php linkTo("command","annotatefunction","annotatefunction");?>(<span class="arg">f</span>, <span class="arg">g</span>, <span class="arg">I</span>, <span class="arg">d</span>, <span class="arg">x0</span>) : (<span class="type">function</span>, <span class="type">function</span>, <span class="type">range</span>, <span class="type">range</span>, <span class="type">constant</span>) -&gt; <span class="type">function</span></span> 
 
</div> 
<div class="divParameters"> 
<h2 class="category">Parameters: </h2> 
<ul> 
<li><span class="arg">f</span> is a function.</li> 
<li><span class="arg">g</span> is a function, in most cases a polynomial.</li> 
<li><span class="arg">I</span> is an interval.</li> 
<li><span class="arg">d</span> is an interval.</li> 
<li><span class="arg">x0</span> is a constant (default value is 0 when not provided).</li> 
</ul> 
</div> 
<div class="divDescription"> 
<h2 class="category">Description: </h2><ul> 
<li>When a given function <span class="arg">f</span> is to be evaluated at several points of a given 
interval <span class="arg">I</span> to a given precision, it might be useful to precompute a good 
approximant <span class="arg">g</span> of <span class="arg">f</span> and further evaluate it instead of <span class="arg">f</span> when the 
approximation is good enough to provide the desire precision. If <span class="arg">f</span> is a 
complicated expression, whereas <span class="arg">g</span> is, e.g., a polynomial of low degree, 
the cost of precomputing <span class="arg">g</span> can be well compensated by the gain of time in 
each subsequent evaluation. The purpose of <?php linkTo("command","annotatefunction","annotatefunction");?> is to provide 
such a mechanism to the user. 
</li><li>When using <?php linkTo("command","annotatefunction","annotatefunction");?>(<span class="arg">f</span>, <span class="arg">g</span>, <span class="arg">I</span>, <span class="arg">d</span>, <span class="arg">x0</span>), 
resp. <?php linkTo("command","annotatefunction","annotatefunction");?>(<span class="arg">f</span>, <span class="arg">g</span>, <span class="arg">I</span>, <span class="arg">d</span>) (where <span class="arg">x0</span> is assumed to be 
zero), it is assumed that 
  
                forall x in I, f(x) - g(x - x0) in d. 
  
It is the user responsibility to ensure this property. Otherwise, any 
subsequent use of <span class="arg">f</span> on points of <span class="arg">I</span> might lead to incorrect values. 
</li><li>A call to <?php linkTo("command","annotatefunction","annotatefunction");?>(<span class="arg">f</span>, <span class="arg">g</span>, <span class="arg">I</span>, <span class="arg">d</span>, <span class="arg">x0</span>) annotates the given 
Sollya function object <span class="arg">f</span> with the approximation <span class="arg">g</span>. In further use, when 
asked to evaluate <span class="arg">f</span> on a point x of <span class="arg">I</span>, Sollya will first evaluate <span class="arg">g</span> 
on x-x0 and check if the result is accurate enough in the given context 
(accounting for the fact that the error of approximation between the true 
value and g(x-x0) belongs to <span class="arg">d</span>). If not (and only in this case), an 
evaluation of the expression of <span class="arg">f</span> on x is performed. 
</li><li>The approximation <span class="arg">g</span> can be any Sollya function but particular 
performance is expected when <span class="arg">g</span> is a polynomial. Upon annotation with a 
polynomial, precomputations are performed to analyze certain properties of 
the given approximation polynomial. 
</li><li><?php linkTo("command","annotatefunction","annotatefunction");?> updates the internal representation of <span class="arg">f</span> so as to 
persistently keep this information attached with the Sollya object 
representing <span class="arg">f</span>. In particular, the annotation is persistent through copy 
or use of <span class="arg">f</span> as a subexpression to build up bigger expressions. Notice 
however, that there is no way of deducing an annotation for the derivative 
of <span class="arg">f</span> from an annotation of <span class="arg">f</span>. So, in general, it should not be expected 
that <?php linkTo("command","diff","diff");?>(<span class="arg">f</span>) will be automatically annotated (notice, however that <span class="arg">f</span> 
might be a subexpression of its derivative, e.g., for <span class="arg">f</span>=<?php linkTo("command","exp","exp");?> or <span class="arg">f</span>=<?php linkTo("command","tan","tan");?>, in 
which case the corresponding subexpressions of the derivative could inherit 
the annotations from <span class="arg">f</span>. It is currently not specified whether Sollya does 
this automatically or not). 
</li><li><?php linkTo("command","annotatefunction","annotatefunction");?> really is an imperative statement that modifies the 
internal representation of <span class="arg">f</span>. However, for convenience <?php linkTo("command","annotatefunction","annotatefunction");?> 
returns <span class="arg">f</span> itself. 
</li><li>Sollya function objects can be annotated more than once with different 
approximations on different domains, that do not need to be disjoint. Upon 
evaluation of the annotated function object, Sollya chooses an 
approximation annotation (if any) that provides for sufficient accuracy at 
the evaluation point. It is not specified in which order Sollya tries 
different possible annotations when several are available for a given 
point <span class="arg">x</span>. 
</ul> 
</div> 
<div class="divExamples"> 
<div class="divExample"> 
<h2 class="category">Example 1: </h2> 
&nbsp;&nbsp;&nbsp;&gt; verbosity=1!;<br> 
&nbsp;&nbsp;&nbsp;&gt; procedure EXP(X,n,p) {<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;var res, oldPrec;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;oldPrec = prec;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;prec = p!;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"Using procedure function exponential with X=" @ X @ ", n=" @ n @ ", and p=" @ p;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;res = exp(X);<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;prec = oldPrec!;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;return res;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;};<br> 
&nbsp;&nbsp;&nbsp;&gt; g = function(EXP);<br> 
&nbsp;&nbsp;&nbsp;&gt; p = 46768052394588893382516870161332864698044514954899b-165 + x * (23384026197294446691258465802074096632225783601255b-164 + x * (5846006549323611672948426613035653821819225877423b-163 + x * (3897337699549074448627696490806815137319821946501b-164 + x * (7794675399098148717422744621371434831048848817417b-167 + x * (24942961277114075921122941174178849425809856036737b-171 + x * (8314320425704876115613838900105097456456371179471b-172 + x * (19004160973039701371579356991645932289422670402995b-176 + x * (19004160972669324148912122254449912156003926801563b-179 + x * (33785175062542597526738679493857229456702396042255b-183 + x * (6757035113643674378393625988264926886191860669891b-184 + x * (9828414707511252769908089206114262766633532289937b-188 + x * (26208861108003813314724515233584738706961162212965b-193 + x * (32257064253325954315953742396999456577223350602741b-197 + x * (578429089657689569703509185903214676926704485495b-195 + x * 2467888542176675658523627105540996778984959471957b-201))))))))))))));<br> 
&nbsp;&nbsp;&nbsp;&gt; h = annotatefunction(g, p, [-1/2;1/2], [-475294848522543b-124;475294848522543b-124]);<br> 
&nbsp;&nbsp;&nbsp;&gt; h == g;<br> 
&nbsp;&nbsp;&nbsp;true<br> 
&nbsp;&nbsp;&nbsp;&gt; prec = 24;<br> 
&nbsp;&nbsp;&nbsp;The precision has been set to 24 bits.<br> 
&nbsp;&nbsp;&nbsp;&gt; h(0.25);<br> 
&nbsp;&nbsp;&nbsp;Warning: rounding has happened. The value displayed is a faithful rounding to 24 bits of the true result.<br> 
&nbsp;&nbsp;&nbsp;1.2840254<br> 
&nbsp;&nbsp;&nbsp;&gt; prec = 165;<br> 
&nbsp;&nbsp;&nbsp;The precision has been set to 165 bits.<br> 
&nbsp;&nbsp;&nbsp;&gt; h(0.25);<br> 
&nbsp;&nbsp;&nbsp;Using procedure function exponential with X=[0.25;0.25], n=0, and p=185<br> 
&nbsp;&nbsp;&nbsp;Warning: rounding has happened. The value displayed is a faithful rounding to 165 bits of the true result.<br> 
&nbsp;&nbsp;&nbsp;1.28402541668774148407342056806243645833628086528147<br> 
</div> 
</div> 
<div class="divSeeAlso"> 
<span class="category">See also: </span><?php linkTo("command","chebyshevform","chebyshevform");?>, <?php linkTo("command","taylorform","taylorform");?>, <?php linkTo("command","remez","remez");?>, <?php linkTo("command","supnorm","supnorm");?>, <?php linkTo("command","infnorm","infnorm");?> 
</div> 
