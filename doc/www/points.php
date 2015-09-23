<a name="points"></a> 
<div class="divName"> 
<h2 class="name">Name:</h2> <?php linkTo("command","points","points");?> 
<span class="smallDescription">controls the number of points chosen by Sollya in certain commands. 
</span> 
</div> 
<div class="divLibraryName"> 
<h2 class="libraryname">Library names:</h2> 
<span class="commandline type">void sollya_lib_set_points_and_print(sollya_obj_t)</span> 
<span class="commandline type">void sollya_lib_set_points(sollya_obj_t)</span> 
<span class="commandline type">sollya_obj_t sollya_lib_get_points()</span> 
</div> 
<div class="divUsage"> 
<h2 class="category">Usage: </h2> 
<span class="commandline"><?php linkTo("command","points","points");?> = <span class="arg">n</span> : <span class="type">integer</span> -&gt; <span class="type">void</span></span> 
<span class="commandline"><?php linkTo("command","points","points");?> = <span class="arg">n</span> ! : <span class="type">integer</span> -&gt; <span class="type">void</span></span> 
<span class="commandline"><?php linkTo("command","points","points");?> : <span class="type">constant</span></span> 
 
</div> 
<div class="divParameters"> 
<h2 class="category">Parameters: </h2> 
<ul> 
<li><span class="arg">n</span> represents the number of points</li> 
</ul> 
</div> 
<div class="divDescription"> 
<h2 class="category">Description: </h2><ul> 
<li><?php linkTo("command","points","points");?> is a global variable. Its value represents the number of points 
used in numerical algorithms of Sollya (namely <?php linkTo("command","dirtyinfnorm","dirtyinfnorm");?>, 
<?php linkTo("command","dirtyintegral","dirtyintegral");?>, <?php linkTo("command","dirtyfindzeros","dirtyfindzeros");?>, <?php linkTo("command","plot","plot");?>). 
</ul> 
</div> 
<div class="divExamples"> 
<div class="divExample"> 
<h2 class="category">Example 1: </h2> 
&nbsp;&nbsp;&nbsp;&gt; f=x^2*sin(1/x);<br> 
&nbsp;&nbsp;&nbsp;&gt; points=10;<br> 
&nbsp;&nbsp;&nbsp;The number of points has been set to 10.<br> 
&nbsp;&nbsp;&nbsp;&gt; dirtyfindzeros(f, [0;1]);<br> 
&nbsp;&nbsp;&nbsp;[|0, 0.31830988618379067153776752674502872406891929148092|]<br> 
&nbsp;&nbsp;&nbsp;&gt; points=100;<br> 
&nbsp;&nbsp;&nbsp;The number of points has been set to 100.<br> 
&nbsp;&nbsp;&nbsp;&gt; dirtyfindzeros(f, [0;1]);<br> 
&nbsp;&nbsp;&nbsp;[|0, 2.4485375860291590118289809749617594159147637806224e-2, 3.9788735772973833942220940843128590508614911435115e-2, 4.5472840883398667362538218106432674866988470211559e-2, 5.3051647697298445256294587790838120678153215246819e-2, 6.3661977236758134307553505349005744813783858296184e-2, 7.957747154594766788444188168625718101722982287023e-2, 0.106103295394596890512589175581676241356306430493638, 0.15915494309189533576888376337251436203445964574046, 0.31830988618379067153776752674502872406891929148092|]<br> 
</div> 
</div> 
<div class="divSeeAlso"> 
<span class="category">See also: </span><?php linkTo("command","dirtyinfnorm","dirtyinfnorm");?>, <?php linkTo("command","dirtyintegral","dirtyintegral");?>, <?php linkTo("command","dirtyfindzeros","dirtyfindzeros");?>, <?php linkTo("command","plot","plot");?>, <?php linkTo("command","diam","diam");?>, <?php linkTo("command","prec","prec");?> 
</div> 
