<a name="annotatefunction"></a> 
<div class="divName"> 
<h2 class="name">Name:</h2> <?php linkTo("command","annotatefunction","annotatefunction");?> 
<span class="smallDescription">Annotates a Sollya function object with an approximation that is faster to evaluate 
</span> 
</div> 
<div class="divLibraryName"> 
<h2 class="libraryname">Library name:</h2> 
<span class="commandline type">sollya_obj_t sollya_lib_annotatefunction(sollya_obj_t, sollya_obj_t, </span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sollya_obj_t, sollya_obj_t, ...);</span> 
</div> 
<div class="divUsage"> 
<h2 class="category">Usage: </h2> 
<span class="commandline"><?php linkTo("command","annotatefunction","annotatefunction");?>(<span class="arg">f</span>, <span class="arg">g</span>, <span class="arg">I</span>, <span class="arg">d</span>) : (<span class="type">function</span>, <span class="type">function</span>, <span class="type">range</span>, <span class="type">range</span>) -&gt; <span class="type">function</span></span> 
<span class="commandline"><?php linkTo("command","annotatefunction","annotatefunction");?>(<span class="arg">f</span>, <span class="arg">g</span>, <span class="arg">I</span>, <span class="arg">d</span>, <span class="arg">t</span>) : (<span class="type">function</span>, <span class="type">function</span>, <span class="type">range</span>, <span class="type">range</span>, <span class="type">constant</span>) -&gt; <span class="type">function</span></span> 
 
</div> 
<div class="divParameters"> 
<h2 class="category">Parameters: </h2> 
<ul> 
<li><span class="arg">f</span> is a function.</li> 
<li><span class="arg">g</span> is a function, in most cases a polynomial.</li> 
<li><span class="arg">I</span> is an interval.</li> 
<li><span class="arg">d</span> is an interval.</li> 
<li><span class="arg">t</span> is a constant.</li> 
</ul> 
</div> 
<div class="divDescription"> 
<h2 class="category">Description: </h2><ul> 
<li><?php linkTo("command","annotatefunction","annotatefunction");?>(<span class="arg">f</span>, <span class="arg">g</span>, <span class="arg">I</span>, <span class="arg">d</span>, <span class="arg">t</span>) resp. <?php linkTo("command","annotatefunction","annotatefunction");?>(<span class="arg">f</span>, <span class="arg">g</span>, <span class="arg">I</span>, <span class="arg">d</span>)  
(where <span class="arg">t</span> is assumed to be zero) annotates the given Sollya 
function object <span class="arg">f</span> with an approximation <span class="arg">g</span> in such a way that upon 
evaluation of <span class="arg">f</span> at a point in <span class="arg">I</span>, <span class="arg">g</span> is used instead of <span class="arg">f</span> as 
long as the error between <span class="arg">f</span> and <span class="arg">g</span>, which must be bounded over all 
<span class="arg">I</span> by <span class="arg">d</span>, is not too large. The approximation <span class="arg">g</span> is not evaluated 
at the same point as <span class="arg">f</span> would be evaluated but after subtraction of  
<span class="arg">t</span> from that point. 
</li><li>In order not to provoke any incorrect behavior, the given <span class="arg">f</span>, <span class="arg">g</span>, 
<span class="arg">I</span>, <span class="arg">d</span> and <span class="arg">t</span> must satisfy the following property: 
forall x in I, f(x) - g(x - t) in d. 
</li><li>The approximation <span class="arg">g</span> can be any Sollya function but particular 
performance is expected when <span class="arg">g</span> is a polynomial. Upon annotation with 
a polynomial, precomputations are performed to analyze certain 
properties of the given approximation polynomial. 
</li><li>The <?php linkTo("command","annotatefunction","annotatefunction");?> command returns the annotated Sollya function 
object. This object is indistinguishable from the original <span class="arg">f</span> with 
respect to Sollya comparisons, printing and so on. Unless Sollya has 
been compiled in some particular debug mode, the <?php linkTo("command","annotatefunction","annotatefunction");?> command does 
not deeply copy the function object <span class="arg">f</span>, so the original function 
object <span class="arg">f</span> also bears the annotation. 
</li><li>Sollya function objects can be annotated more than once with 
different approximations on different domains, that do not need to be 
disjoint. Upon evaluation of the annotated function object, Sollya 
chooses the approximation annotation that provides for sufficient 
accuracy at the evaluation point. 
</li><li>As an absolute error bound <span class="arg">d</span> is given for annotation, Sollya cannot 
reuse an approximation annotation for the derivative of the annotated 
function. If any (higher) derivative of a function needs to bear an 
annotation, that derivative needs to be annotated using <?php linkTo("command","annotatefunction","annotatefunction");?>. 
</ul> 
</div> 
<div class="divExamples"> 
<div class="divExample"> 
<h2 class="category">Example 1: </h2> 
&nbsp;&nbsp;&nbsp;&gt; procedure EXP(X,n,p) {<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;var res, oldPrec;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;oldPrec = prec;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;prec = p!;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"Using procedure function exponential";<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;res = exp(X);<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;prec = oldPrec!;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;return res;<br> 
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;};<br> 
&nbsp;&nbsp;&nbsp;&gt; g = function(EXP);<br> 
&nbsp;&nbsp;&nbsp;&gt; p = 46768052394588893382516870161332864698044514954899b-165 + x * (23384026197294446691258465802074096632225783601255b-164 + x * (5846006549323611672948426613035653821819225877423b-163 + x * (3897337699549074448627696490806815137319821946501b-164 + x * (7794675399098148717422744621371434831048848817417b-167 + x * (24942961277114075921122941174178849425809856036737b-171 + x * (8314320425704876115613838900105097456456371179471b-172 + x * (19004160973039701371579356991645932289422670402995b-176 + x * (19004160972669324148912122254449912156003926801563b-179 + x * (33785175062542597526738679493857229456702396042255b-183 + x * (6757035113643674378393625988264926886191860669891b-184 + x * (9828414707511252769908089206114262766633532289937b-188 + x * (26208861108003813314724515233584738706961162212965b-193 + x * (32257064253325954315953742396999456577223350602741b-197 + x * (578429089657689569703509185903214676926704485495b-195 + x * 2467888542176675658523627105540996778984959471957b-201))))))))))))));<br> 
&nbsp;&nbsp;&nbsp;&gt; h = annotatefunction(g, p, [-1/2;1/2], [-475294848522543b-124;475294848522543b-124]);<br> 
&nbsp;&nbsp;&nbsp;&gt; prec = 24;<br> 
&nbsp;&nbsp;&nbsp;The precision has been set to 24 bits.<br> 
&nbsp;&nbsp;&nbsp;&gt; "h(0.25) = ";<br> 
&nbsp;&nbsp;&nbsp;h(0.25) = <br> 
&nbsp;&nbsp;&nbsp;&gt; h(0.25);<br> 
&nbsp;&nbsp;&nbsp;1.2840254<br> 
&nbsp;&nbsp;&nbsp;&gt; prec = 72;<br> 
&nbsp;&nbsp;&nbsp;The precision has been set to 72 bits.<br> 
&nbsp;&nbsp;&nbsp;&gt; "h(0.25) = ";<br> 
&nbsp;&nbsp;&nbsp;h(0.25) = <br> 
&nbsp;&nbsp;&nbsp;&gt; h(0.25);<br> 
&nbsp;&nbsp;&nbsp;Using procedure function exponential<br> 
&nbsp;&nbsp;&nbsp;Using procedure function exponential<br> 
&nbsp;&nbsp;&nbsp;1.2840254166877414840735<br> 
</div> 
</div> 
<div class="divSeeAlso"> 
<span class="category">See also: </span><?php linkTo("command","taylorform","taylorform");?>, <?php linkTo("command","remez","remez");?>, <?php linkTo("command","supnorm","supnorm");?>, <?php linkTo("command","infnorm","infnorm");?> 
</div> 
