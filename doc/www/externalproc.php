<a name="externalproc"></a> 
<div class="divName"> 
<h2 class="name">Name:</h2> <?php linkTo("command","externalproc","externalproc");?> 
<span class="smallDescription">binds an external code to a Sollya procedure 
</span> 
</div> 
<div class="divLibraryName"> 
<h2 class="libraryname">Library names:</h2> 
<span class="commandline type">&nbsp;&nbsp;sollya_obj_t sollya_lib_externalprocedure(sollya_externalprocedure_type_t, </span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sollya_externalprocedure_type_t *,</span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;int, char *, void *);</span> 
<span class="commandline type">&nbsp;&nbsp;sollya_obj_t sollya_lib_externalprocedure_with_data(</span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sollya_externalprocedure_type_t, </span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sollya_externalprocedure_type_t *,</span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;int, char *, void *, void *,</span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;void (*)(void *));</span> 
</div> 
<div class="divUsage"> 
<h2 class="category">Usage: </h2> 
<span class="commandline"><?php linkTo("command","externalproc","externalproc");?>(<span class="arg">identifier</span>, <span class="arg">filename</span>, <span class="arg">argumenttype</span> -&gt; <span class="arg">resulttype</span>) : (<span class="type">identifier type</span>, <span class="type">string</span>, <span class="type">type type</span>, <span class="type">type type</span>) -&gt; <span class="type">void</span></span> 
 
</div> 
<div class="divParameters"> 
<h2 class="category">Parameters: </h2> 
<ul> 
<li><span class="arg">identifier</span> represents the identifier the code is to be bound to</li> 
<li><span class="arg">filename</span> of type <span class="type">string</span> represents the name of the object file where the code of procedure can be found</li> 
<li><span class="arg">argumenttype</span> represents a definition of the types of the arguments of the Sollya procedure and the external code</li> 
<li><span class="arg">resulttype</span> represents a definition of the result type of the external code</li> 
</ul> 
</div> 
<div class="divDescription"> 
<h2 class="category">Description: </h2><ul> 
<li><?php linkTo("command","externalproc","externalproc");?> allows for binding the Sollya identifier <span class="arg">identifier</span> to an 
external code. After this binding, when Sollya encounters <span class="arg">identifier</span> 
applied to a list of actual parameters, it will evaluate these parameters and 
call the external code with these parameters. If the external code indicated 
success, it will receive the result produced by the external code, transform 
it to Sollya's internal representation and return it. 
<br><br> 
In order to allow correct evaluation and typing of the data in parameter and 
in result to be passed to and received from the external code, <?php linkTo("command","externalproc","externalproc");?> 
has a third parameter <span class="arg">argumenttype</span> -&gt; <span class="arg">resulttype</span>. Both <span class="arg">argumenttype</span> and 
<span class="arg">resulttype</span> are one of <?php linkTo("command","void","void");?>, <?php linkTo("command","constant","constant");?>, <?php linkTo("command","function","function");?>, <?php linkTo("command","object","object");?>, <?php linkTo("command","range","range");?>, <?php linkTo("command","integer","integer");?>, 
<?php linkTo("command","string","string");?>, <?php linkTo("command","boolean","boolean");?>, <?php linkTo("command","listof","list of");?> <?php linkTo("command","constant","constant");?>, <?php linkTo("command","listof","list of");?> <?php linkTo("command","function","function");?>, <?php linkTo("command","listof","list of");?> <?php linkTo("command","object","object");?>, 
<?php linkTo("command","listof","list of");?> <?php linkTo("command","range","range");?>, <?php linkTo("command","listof","list of");?> <?php linkTo("command","integer","integer");?>, <?php linkTo("command","listof","list of");?> <?php linkTo("command","string","string");?>, <?php linkTo("command","listof","list of");?> <?php linkTo("command","boolean","boolean");?>. 
<br><br> 
It is worth mentionning that the difference between the data and 
result type <?php linkTo("command","function","function");?> and the type <?php linkTo("command","object","object");?> is minimal and due to 
support of legacy Sollya code. Both Sollya functions and Sollya 
objects are transferred from and to the external procedure thru the C 
type sollya_obj_t. The difference is that 
Sollya will check that a certain object is a mathematical function 
when <?php linkTo("command","function","function");?> is used as a type, and will skip this test if the 
<?php linkTo("command","object","object");?> type is used. Similarly, Sollya relies on an object produced 
by the external procedure to be a mathematical function when <?php linkTo("command","function","function");?> 
is used and will not make this assumption for <?php linkTo("command","object","object");?>. 
<br><br> 
If upon a usage of a procedure bound to an external procedure the type of the 
actual parameters given or its number is not correct, Sollya produces a type 
error. An external function not applied to arguments represents itself and 
prints out with its argument and result types. 
<br><br> 
The external function is supposed to return an integer indicating success. It 
returns its result depending on its Sollya result type as follows. Here, the 
external procedure is assumed to be implemented as a C function.<ul> 
  <li> If the Sollya result type is void, the C function has no pointer 
     argument for the result. 
  </li><li> If the Sollya result type is <?php linkTo("command","constant","constant");?>, the first argument of the 
     C function is of C type mpfr_t *, the result is returned by affecting 
     the MPFR variable. 
  </li><li> If the Sollya result type is <?php linkTo("command","function","function");?>, the first argument of the 
     C function is of C type sollya_obj_t *, the result is returned by 
     affecting the sollya_obj_t variable. 
  </li><li> If the Sollya result type is <?php linkTo("command","object","object");?>, the first argument of the 
     C function is of C type sollya_obj_t *, the result is returned by 
     affecting the sollya_obj_t variable. 
  </li><li> If the Sollya result type is <?php linkTo("command","range","range");?>, the first argument of the C function 
     is of C type mpfi_t *, the result is returned by affecting the MPFI 
     variable. 
  </li><li> If the Sollya result type is <?php linkTo("command","integer","integer");?>, the first argument of the 
     C function is of C type int *, the result is returned by affecting the 
     int variable. 
  </li><li> If the Sollya result type is <?php linkTo("command","string","string");?>, the first argument of the 
     C function is of C type char **, the result is returned by the char * 
     pointed with a new char *. 
  </li><li> If the Sollya result type is <?php linkTo("command","boolean","boolean");?>, the first argument of the 
     C function is of C type int *, the result is returned by affecting the 
     int variable with a boolean value. 
  </li><li> If the Sollya result type is <?php linkTo("command","listof","list of");?> type, the first argument of the 
     C function is of a C type depending on the Sollya return type:<ul> 
       <li> For a list of <?php linkTo("command","constant","constant");?>: sollya_constant_list_t * 
       </li><li> For a list of <?php linkTo("command","function","function");?>: sollya_obj_list_t * 
       </li><li> For a list of <?php linkTo("command","object","object");?>: sollya_obj_list_t * 
       </li><li> For a list of <?php linkTo("command","range","range");?>: sollya_constant_list_t * 
       </li><li> For a list of <?php linkTo("command","integer","integer");?>: sollya_int_list_t * 
       </li><li> For a list of <?php linkTo("command","string","string");?>: sollya_string_list_t * 
       </li><li> For a list of <?php linkTo("command","boolean","boolean");?>: sollya_boolean_list_t * </li></ul> 
</li></ul> 
The external procedure affects its possible pointer argument if and only if 
it succeeds. This means, if the function returns an integer indicating 
failure, it does not leak any memory to the encompassing environment. 
<br><br> 
The external procedure receives its arguments as follows: If the Sollya 
argument type is <?php linkTo("command","void","void");?>, no argument array is given. Otherwise the C function 
receives a C void ** argument representing an array of size equal to the 
arity of the function where each entry (of C type void *) represents a value 
with a C type depending on the corresponding Sollya type.<ul> 
  <li> If the Sollya type is <?php linkTo("command","constant","constant");?>, the void * is to be cast to mpfr_t *. 
  </li><li> If the Sollya type is <?php linkTo("command","function","function");?>, the void * is to be cast to 
     sollya_obj_t. 
  </li><li> If the Sollya type is <?php linkTo("command","object","object");?>, the void * is to be cast to sollya_obj_t. 
  </li><li> If the Sollya type is <?php linkTo("command","range","range");?>, the void * is to be cast to mpfi_t *. 
  </li><li> If the Sollya type is <?php linkTo("command","integer","integer");?>, the void * is to be cast to int *. 
  </li><li> If the Sollya type is <?php linkTo("command","string","string");?>, the void * is to be cast to char *. 
  </li><li> If the Sollya type is <?php linkTo("command","boolean","boolean");?>, the void * is to be cast to int *. 
  </li><li> If the Sollya type is <?php linkTo("command","listof","list of");?> type, the void * is to be cast to a list 
     of a type depending on the type of the list argument:<ul> 
       <li> For a list of <?php linkTo("command","constant","constant");?>: sollya_constant_list_t 
       </li><li> For a list of <?php linkTo("command","function","function");?>: sollya_obj_list_t 
       </li><li> For a list of <?php linkTo("command","object","object");?>: sollya_obj_list_t 
       </li><li> For a list of <?php linkTo("command","range","range");?>: sollya_interval_list_t 
       </li><li> For a list of <?php linkTo("command","integer","integer");?>: sollya_int_list_t 
       </li><li> For a list of <?php linkTo("command","string","string");?>: sollya_string_list_t 
       </li><li> For a list of <?php linkTo("command","boolean","boolean");?>: sollya_boolean_list_t </li></ul> 
</li></ul> 
The external procedure is not supposed to alter the memory pointed by its 
array argument void **. 
<br><br> 
In both directions (argument and result values), empty lists are represented 
by NULL pointers. 
<br><br> 
Similarly to internal procedures, externally bounded procedures can be 
considered to be objects inside Sollya that can be assigned to other 
variables, stored in list etc. 
</li><li>The user should be aware that they may use the Sollya library in external 
codes to be dynamically bound to Sollya using <?php linkTo("command","externalproc","externalproc");?>. On most systems, 
it suffices to include the header of the Sollya library into the source code 
of the external procedure. Linking with the actual Sollya library is not 
necessary on most systems; as the interactive Sollya executable contains a 
superset of the Sollya library functions. On some systems, linking with the 
Sollya library or some of its dependencies may be necessary. 
<br><br> 
In particular, the Sollya library &ndash;&nbsp;and, of course, its header file&nbsp;&ndash; 
contain a certain set of functions to manipulate lists with elements of 
certain types, such as sollya_constant_list_t, sollya_obj_list_t and so on. 
As explained above, these types are passed in argument to (and received back 
thru a reference from) an external procedure. These list manipulation 
functions are not strictly necessary to the use of the Sollya library in 
free-standing applications that do not use the functionality provided with 
<?php linkTo("command","externalproc","externalproc");?>. They are therefore provided as-is without any further 
documentation, besides the comments given in the Sollya library header file. 
</li><li>The dynamic object file whose name is given to <?php linkTo("command","externalproc","externalproc");?> for binding of 
an external procedure may also define a destructor function 
int sollya_external_lib_close(void). If Sollya finds such a destructor 
function in the dynamic object file, it will call that function when closing 
the dynamic object file again. This happens when Sollya is terminated or when 
the current Sollya session is restarted using <?php linkTo("command","restart","restart");?>. The purpose of the 
destructor function is to allow the dynamically bound code to free any memory 
that it might have allocated before Sollya is terminated or restarted. 
<br><br> 
The dynamic object file is not necessarily needed to define a destructor 
function. This ensure backward compatibility with older Sollya external 
library function object files. 
<br><br> 
When defined, the destructor function is supposed to return an integer 
value indicating if an error has happened. Upon success, the destructor 
functions is to return a zero value, upon error a non-zero value. 
</ul> 
</div> 
<div class="divExamples"> 
<div class="divExample"> 
<h2 class="category">Example 1: </h2> 
&nbsp;&nbsp;&nbsp;&gt; bashexecute("gcc -fPIC -Wall -c externalprocexample.c");<br> 
&nbsp;&nbsp;&nbsp;&gt; bashexecute("gcc -fPIC -shared -o externalprocexample externalprocexample.o");<br> 
&nbsp;&nbsp;&nbsp;&gt; externalproc(foo, "./externalprocexample", (integer, integer) -&gt; integer);<br> 
&nbsp;&nbsp;&nbsp;&gt; foo;<br> 
&nbsp;&nbsp;&nbsp;foo<br> 
&nbsp;&nbsp;&nbsp;&gt; foo(5, 6);<br> 
&nbsp;&nbsp;&nbsp;11<br> 
&nbsp;&nbsp;&nbsp;&gt; verbosity = 1!;<br> 
&nbsp;&nbsp;&nbsp;&gt; foo();<br> 
&nbsp;&nbsp;&nbsp;Warning: at least one of the given expressions or a subexpression is not correctly typed<br> 
&nbsp;&nbsp;&nbsp;or its evaluation has failed because of some error on a side-effect.<br> 
&nbsp;&nbsp;&nbsp;error<br> 
&nbsp;&nbsp;&nbsp;&gt; a = foo;<br> 
&nbsp;&nbsp;&nbsp;&gt; a(5,6);<br> 
&nbsp;&nbsp;&nbsp;11<br> 
</div> 
</div> 
<div class="divSeeAlso"> 
<span class="category">See also: </span><?php linkTo("command","library","library");?>, <?php linkTo("command","libraryconstant","libraryconstant");?>, <?php linkTo("command","externalplot","externalplot");?>, <?php linkTo("command","bashexecute","bashexecute");?>, <?php linkTo("command","void","void");?>, <?php linkTo("command","constant","constant");?>, <?php linkTo("command","function","function");?>, <?php linkTo("command","range","range");?>, <?php linkTo("command","integer","integer");?>, <?php linkTo("command","string","string");?>, <?php linkTo("command","boolean","boolean");?>, <?php linkTo("command","listof","list of");?> 
</div> 
