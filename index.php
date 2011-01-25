<html>
  <head>
    <title>Cello Unit Tests</title>
    <link href="www/cello.css" rel="stylesheet" type="text/css">
  </head>
<body>
  <h1>Cello Unit Tests</h1>

<ul>
  <li><a href="http://lca.ucsd.edu/projects/cello/">Cello User Site</a></li>
  <li><a href="http://client65-77.sdsc.edu/cello/">Cello Developer Site</a></li>

</ul>
<?php

function component($component) {
   printf ("<hr><h2>$component Component</h2><a name=\"$component\">");
}

function tests($component,$testrun,$output) {

   $parallel_types = array("serial","mpi","ampi","charm");

   $source_file = "src/$component/$testrun.cpp";
   $source_html = "<a href=\"$source_file\">$testrun.cpp</a>";

   echo "<h3>Test: $source_html</h3>\n";
   if (! file_exists($source_file)) {
     echo "<p><strong class=fail>$source_file does not exist</strong></p>";
   } else {
     echo "<table>\n";
     echo "<tr>\n";
     echo "<th></th>";
     for ($i = 0; $i<sizeof($parallel_types); ++$i) {
       echo "<th>$parallel_types[$i]</th>";
     }
     echo "</tr>\n";

     echo "<tr>\n";
     //--------------------------------------------------
     echo "<th>Output</th>";
     //--------------------------------------------------
     for ($i = 0; $i<sizeof($parallel_types); ++$i) {
       $output_file = "test/$parallel_types[$i]-$output.unit";
       if (file_exists($output_file)) {
	 $output_html = "<a href=\"$output_file\">$output.unit</a>";
	 system("awk 'BEGIN{c=0}; /UNIT TEST END/ {c=1}; END{ if (c==0) print \"<td class=fail><a href='$output_file'>incomplete</a></td>\"; if (c!=0) print \"<td class=pass><a href='$output_file'>complete</a></td>\"}' < $output_file");
       } else {
	 echo "<td></td>";
       }
     }
     echo "</tr>\n";

     echo "<tr>\n";
     //--------------------------------------------------
     echo "<th>Date</th>";
     //--------------------------------------------------
     for ($i = 0; $i<sizeof($parallel_types); ++$i) {
       $output_file = "test/$parallel_types[$i]-$output.unit";
       if (file_exists($output_file)) {
	 $output_html = date ("Y-m-d", filemtime($output_file));
	 echo "<td class=pass>$output_html</td>";
       } else {
	 echo "<td class=pass></td>";
       }
     }
     echo "</tr>\n";

     echo "<tr>\n";
     //--------------------------------------------------
     echo "<th>Time</th>";
     //--------------------------------------------------
     for ($i = 0; $i<sizeof($parallel_types); ++$i) {
       $output_file = "test/$parallel_types[$i]-$output.unit";
       if (file_exists($output_file)) {
	 $output_html = date ("H:i:s", filemtime($output_file));
	 echo "<td class=pass>$output_html</td>";
       } else {
	 echo "<td class=pass></td>";
       }
     }
     echo "</tr>\n";

     echo "<tr>\n";
     //--------------------------------------------------
     echo "<th>Passed</th>";
     //--------------------------------------------------
     for ($i = 0; $i<sizeof($parallel_types); ++$i) {
       $output_file = "test/$parallel_types[$i]-$output.unit";
       if (file_exists($output_file)) {
	 echo "<td class=pass>";
	 system("grep pass $output_file | wc -l");
	 echo "</td>";
       } else {
	 echo "<td class=pass></td>";
       }
     }
     echo "</tr>\n";

     echo "<tr>\n";
     //--------------------------------------------------
     echo "<th>Failed</th>";
     //--------------------------------------------------
     for ($i = 0; $i<sizeof($parallel_types); ++$i) {
       $output_file = "test/$parallel_types[$i]-$output.unit";
       if (file_exists($output_file)) {
	 system("awk 'BEGIN{c=0}; /FAIL/ {c=c+1}; END {if (c!=0) print\"<td class=fail>\",c,\"</td>\"; if (c==0) print \"<td class=pass></td>\"}' < $output_file");
       } else {
	 echo "<td class=pass></td>";
       }
     }
     echo "</tr>\n";

     echo "</tr></table></br/>\n";

     echo "<table><tr>";
     for ($i = 0; $i<sizeof($parallel_types); ++$i) {
       $parallel_type = $parallel_types[$i];
       $output_file = "test/$parallel_type-$output.unit";
       if (file_exists($output_file)) {
	 test($parallel_type,$testrun,"FAIL");
       }
     }
     echo "</tr></table></br/>";
   }
 };

function test($parallel_type,$testrun,$type) {
  $ltype = strtolower($type);
  if ($type == "pass") {
    $cols = "$4,$6";
    $itemtext  = "";
    $rowtext = "</tr><tr>";
  } else {
    $cols = "\"$file\",$4,$6";
    $itemtext  = "";
    $rowtext = "</tr><tr>";
  }

  $output = "test/$parallel_type-$testrun.unit";
  $count = exec("grep $type $output | wc -l");
  if ($count == 0) {
#     echo "<strong >no ${ltype}ed tests</strong></br/>";
  } else {
     echo "<th class=$type><strong>$parallel_type ${ltype}ed</strong></th> ";
     system ("awk 'BEGIN {c=1}; / $type /{split($3,a,\"\/\"); print \"<td class=$type> \",$cols , \" </td>$itemtext\"; c=c+1}; {if (c==5) {c=0; print \"$rowtext\"}}' < $output");
     echo "</tr><tr></tr>";
  }
     
};

?>

<!-- <?php 
   if (file_exists('CELLO_PLATFORM')) {
     echo "<h2>Platform tested: ";
     system ("cat 'CELLO_PLATFORM'");
     echo "</h2>\n";
   } else {
     echo "<h2>Unknown platform</h2>\n";
   }
?> -->

<p>
Test data shown on this page is automatically generated whenever <code>Cello</code> is compiled.  FAIL status usually indicates that the corresponding function has not been implemented yet.  Empty tables indicate that the unit test files have been deleted, probably with "<code>scons -c</code>" (<code>scons</code> version of "<code>make clean</code>").
</p>

<hr>

<h2>Summary</h2>

<?php

function test_summary($component,$test_output)
{
  printf ("<tr><th><a href=\"#$component\">$component</a></th>\n");

  $parallel_types  = array("serial","mpi","ampi","charm");

  for ($i = 0; $i<sizeof($parallel_types); ++$i) {

    $output_files = "";
    for ($test = 0; $test<sizeof($test_output); ++$test) {
      $output = $test_output[$test];
      $output_files = "$output_files test/$parallel_types[$i]-test_$output.unit";
    }
    system("cat $output_files | awk 'BEGIN {c=0}; /pass/{c=c+1}; END{if (c==0) {print \"<td></td>\"} else {print \"<td class=pass>\",c,\"</td>\";}} '");

  }


  for ($i = 0; $i<sizeof($parallel_types); ++$i) {

    $output_files = "";
    for ($test = 0; $test<sizeof($test_output); ++$test) {
      $output = $test_output[$test];
      $output_files = "$output_files test/$parallel_types[$i]-test_$output.unit";
    }

    system("cat $output_files | awk 'BEGIN {c=0}; /FAIL/{c=c+1}; END{if (c==0) {print \"<td></td>\"} else {print \"<td class=fail>\",c,\"</td>\";}} '");
  }

  printf ("</tr>\n");
}
?>


<?php
printf ("<table><tr><th></th>\n");
$parallel_labels = array("serial", "mpi",  "ampi",   "charm");
for ($k = 0; $k < 2; $k ++) {
  for ($i = 0; $i < sizeof($parallel_labels); ++$i) {
    printf ("<th> $parallel_labels[$i] </th>");
  }
 }

test_summary("Disk",array("FileHdf5","FileIfrit")); 
test_summary("Distribute",array("")); 
test_summary("Error",array("Error")); 
test_summary("Enzo",array("enzo-p","ppm_image","ppm_implosion","ppm_implosion3","ppml_blast","ppml_implosion")); 
test_summary("Field",array("FieldBlock","FieldDescr","FieldFaces")); 
test_summary("Memory",array("Memory")); 
test_summary("Mesh",array("DataDescr","Mesh","Patch")); 
test_summary("Method",array("")); 
test_summary("Monitor",array("Monitor")); 
test_summary("Parallel",array("GroupProcessMpi","Layout")); 
test_summary("Parameters",array("Parameters")); 
test_summary("Particles",array("")); 
test_summary("Performance",array("Performance")); 
test_summary("Portal",array("")); 
test_summary("Schedule",array("Schedule")); 
test_summary("Simulation",array("Simulation")); ?>
</tr></table>

</ul>




<?php

component("Disk");

tests("Disk","test_FileHdf5", "test_FileHdf5");
tests("Disk","test_FileIfrit","test_FileIfrit");

component("Distribute");

component("Error");
tests("Error","test_Error","test_Error");


component("Field");
tests("Field","test_FieldDescr","test_FieldDescr");
tests("Field","test_FieldBlock","test_FieldBlock");
tests("Field","test_FieldFaces","test_FieldFaces");

component("Memory");
tests("Memory","test_Memory","test_Memory");


component("Mesh");

tests("Mesh","test_DataDescr","test_DataDescr"); 
tests("Mesh","test_Mesh","test_Mesh"); 
tests("Mesh","test_Patch","test_Patch"); 

component("Method");

component("Monitor");
tests("Monitor","test_Monitor","test_Monitor");
printf ("<img src=\"monitor_image_1.png\"></img>\n");
printf ("<img src=\"monitor_image_2.png\"></img>\n");
printf ("<img src=\"monitor_image_3.png\"></img>\n");
printf ("<img src=\"monitor_image_4.png\"></img>\n");

component("Parallel");
tests("Parallel","test_GroupProcessMpi","test_GroupProcessMpi");
tests("Parallel","test_Layout","test_Layout");


component("Parameters");
tests("Parameters","test_Parameters","test_Parameters");

component("Particles");

component("Performance");
tests("Performance","test_Performance","test_Performance");

component("Portal");

component("Schedule");
tests("Schedule","test_Schedule","test_Schedule");

component("Simulation");
tests("Simulation","test_Simulation","test_Simulation");

component("Task");

component("Test");

component("Enzo");

?>


<h4>enzo-p</h4>

<?php tests("Enzo","enzo-p","enzo-p_implosion"); ?>


<h4>test_ppm ppm_image</h4>

<?php tests("Enzo/ppm","test_ppm","test_ppm_image"); ?>

<table>
<tr><th>cycle = 0</th><th>cycle=10</th><th>cycle=30</th></tr>
<tr>
<td><img width=320 src="slice-ppm-image-000000.png"></img></td>
<td><img width=320  src="slice-ppm-image-000010.png"></img></td>
<td><img width=320  src="slice-ppm-image-000030.png"></img></td>
</tr>
</table>

<h4>test_ppm ppm_implosion</h4>

<?php tests("Enzo/ppm","test_ppm","test_ppm_implosion"); ?>

<table>
<tr><th>cycle = 0</th><th>cycle=10</th></tr>
<tr>
<td><img width=320 src="slice-ppm-implosion-000000.png"></img></td>
<td><img width=320 src="slice-ppm-implosion-000010.png"></img></td>
</tr>
</table>

<h4>test_ppm ppm_implosion3</h4>

<?php tests("Enzo/ppm","test_ppm","test_ppm_implosion3"); ?>

<table>
<tr>
<th>cycle = 0 project = X</th>
<th>cycle = 0 project = Y</th>
<th>cycle = 0 project = Z</th>
<th>cycle = 10 project = X</th>
<th>cycle = 10 project = Y</th>
<th>cycle = 10 project = Z</th>
</tr>
<tr>
<td><img width=200 src="project-ppm-implosion3-000000-x.png"></img></td>
<td><img width=200 src="project-ppm-implosion3-000000-y.png"></img></td>
<td><img width=200 src="project-ppm-implosion3-000000-z.png"></img></td>
<td><img width=200 src="project-ppm-implosion3-000010-x.png"></img></td>
<td><img width=200 src="project-ppm-implosion3-000010-y.png"></img></td>
<td><img width=200 src="project-ppm-implosion3-000010-z.png"></img></td>
</tr>
</table>

<h3>test_ppml ppml_blast</h3>

<?php tests("Enzo/ppml","test_ppml","test_ppml_blast"); ?>

<table>
<tr>
<th>cycle = 0 project = X</th>
<th>cycle = 0 project = Y</th>
<th>cycle = 0 project = Z</th>
<th>cycle = 10 project = X</th>
<th>cycle = 10 project = Y</th>
<th>cycle = 10 project = Z</th>
</tr>
<tr>
<td><img width=200 src="project-ppml-blast-000000-x.png"></img></td>
<td><img width=200 src="project-ppml-blast-000000-y.png"></img></td>
<td><img width=200 src="project-ppml-blast-000000-z.png"></img></td>
<td><img width=200 src="project-ppml-blast-000010-x.png"></img></td>
<td><img width=200 src="project-ppml-blast-000010-y.png"></img></td>
<td><img width=200 src="project-ppml-blast-000010-z.png"></img></td>
</tr>
</table>

<h3>test_ppml ppml_implosion</h3>

<?php tests("Enzo/ppml","test_ppml","test_ppml_implosion"); ?>

<table>
<tr>
<th>cycle = 0 project = X</th>
<th>cycle = 0 project = Y</th>
<th>cycle = 0 project = Z</th>
<th>cycle = 10 project = X</th>
<th>cycle = 10 project = Y</th>
<th>cycle = 10 project = Z</th>
</tr>
<tr>
<td><img width=200 src="project-ppml-implosion3-000000-x.png"></img></td>
<td><img width=200 src="project-ppml-implosion3-000000-y.png"></img></td>
<td><img width=200 src="project-ppml-implosion3-000000-z.png"></img></td>
<td><img width=200 src="project-ppml-implosion3-000010-x.png"></img></td>
<td><img width=200 src="project-ppml-implosion3-000010-y.png"></img></td>
<td><img width=200 src="project-ppml-implosion3-000010-z.png"></img></td>
</tr>
</table>

<h3>EnzoMethodPpm</h3>

<?php tests("Enzo","test_EnzoMethodPpm","test_EnzoMethodPpm"); ?>

<table>
<tr>
<th></th>
<th>density</th>
<th>velocity_x</th>
<th>velocity_y</th>
<th>total_energy</th></tr>
<tr>
<th>cycle = 0</th>
<td><img width=200 src="ppm-density-00000.png"></img></td>
<td><img width=200 src="ppm-velocity_x-00000.png"></img></td>
<td><img width=200 src="ppm-velocity_y-00000.png"></img></td>
<td><img width=200 src="ppm-total_energy-00000.png"></img></td>
</tr>
<tr>
<th>cycle = 100</th>
<td><img width=200 src="ppm-density-00100.png"></img></td>
<td><img width=200 src="ppm-velocity_x-00100.png"></img></td>
<td><img width=200 src="ppm-velocity_y-00100.png"></img></td>
<td><img width=200 src="ppm-total_energy-00100.png"></img></td>
</tr>
</table>

<h3>TreeK-D2-R2-L?</h3>

<?php tests("Mesh","test_TreeK","test_TreeK-D2-R2-L6"); ?>
<?php tests("Mesh","test_TreeK","test_TreeK-D2-R2-L7"); ?>
<?php tests("Mesh","test_TreeK","test_TreeK-D2-R2-L8"); ?>
<?php tests("Mesh","test_TreeK","test_TreeK-D2-R2-L9"); ?>
<?php tests("Mesh","test_TreeK","test_TreeK-D2-R2-L10"); ?>

<table>
<tr>
<th></th>
<th>levels = 6</th>
<th>levels = 7</th>
<th>levels = 8</th>
<th>levels = 9</th>
<th>levels = 10</th>
</tr>
<tr>
<th>coalesce = false</th>
<td><img width=257 src="TreeK-D=2-R=2-L=6-0.png"></img></td>
<td><img width=257 src="TreeK-D=2-R=2-L=7-0.png"></img></td>
<td><img width=257 src="TreeK-D=2-R=2-L=8-0.png"></img></td>
<td><img width=257 src="TreeK-D=2-R=2-L=9-0.png"></img></td>
<td><img width=257 src="TreeK-D=2-R=2-L=10-0.png"></img></td>
</tr>
<tr>
<th>coalesce = true</th>
<td><img width=257 src="TreeK-D=2-R=2-L=6-1.png"></img></td>
<td><img width=257 src="TreeK-D=2-R=2-L=7-1.png"></img></td>
<td><img width=257 src="TreeK-D=2-R=2-L=8-1.png"></img></td>
<td><img width=257 src="TreeK-D=2-R=2-L=9-1.png"></img></td>
<td><img width=257 src="TreeK-D=2-R=2-L=10-1.png"></img></td>
</tr>
</table>

<h3>TreeK-D2-R4-L?</h3>

<?php tests("Mesh","test_TreeK","test_TreeK-D2-R4-L6"); ?>
<?php tests("Mesh","test_TreeK","test_TreeK-D2-R4-L8"); ?>
<?php tests("Mesh","test_TreeK","test_TreeK-D2-R4-L10"); ?>

<table>
<tr>
<th></th>
<th>levels = 6</th>
<th>levels = 8</th>
<th>levels = 10</th>
</tr>
<tr>
<th>coalesce = false</th>
<td><img width=257 src="TreeK-D=2-R=4-L=6-0.png"></img></td>
<td><img width=257 src="TreeK-D=2-R=4-L=8-0.png"></img></td>
<td><img width=257 src="TreeK-D=2-R=4-L=10-0.png"></img></td>
</tr>
<tr>
<th>coalesce = true</th>
<td><img width=257 src="TreeK-D=2-R=4-L=6-1.png"></img></td>
<td><img width=257 src="TreeK-D=2-R=4-L=8-1.png"></img></td>
<td><img width=257 src="TreeK-D=2-R=4-L=10-1.png"></img></td>
</tr>
</table>

<h3>TreeK-D3-R2-L?</h3>

<?php tests("Mesh","test_TreeK","test_TreeK-D3-R2-L4"); ?>
<?php tests("Mesh","test_TreeK","test_TreeK-D3-R2-L5"); ?>
<?php tests("Mesh","test_TreeK","test_TreeK-D3-R2-L6"); ?>
<?php tests("Mesh","test_TreeK","test_TreeK-D3-R2-L7"); ?>
<?php tests("Mesh","test_TreeK","test_TreeK-D3-R2-L8"); ?>

<table>
<tr>
<th>coalesce = false</th>
<th>levels = 4</th>
<th>levels = 5</th>
<th>levels = 6</th>
<th>levels = 7</th>
<th>levels = 8</th>
</tr>
<tr>
<th>project = X</th>
<td><img width=129 src="TreeK-D=3-R=2-L=4-x-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=5-x-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=6-x-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=7-x-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=8-x-0.png"></img></td>
</tr>
<tr>
<th>project = Y</th>
<td><img width=129 src="TreeK-D=3-R=2-L=4-y-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=5-y-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=6-y-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=7-y-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=8-y-0.png"></img></td>
</tr>
<tr>
<th>project = Z</th>
<td><img width=129 src="TreeK-D=3-R=2-L=4-z-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=5-z-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=6-z-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=7-z-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=8-z-0.png"></img></td>
</tr>
</table>

<p>

<table>
<tr>
<th>coalesce = true</th>
<th>levels = 4</th>
<th>levels = 5</th>
<th>levels = 6</th>
<th>levels = 7</th>
<th>levels = 8</th>
</tr>
<tr>
<th>project = X</th>
<td><img width=129 src="TreeK-D=3-R=2-L=4-x-1.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=5-x-1.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=6-x-1.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=7-x-1.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=8-x-1.png"></img></td>
</tr>
<tr>
<th>project = Y</th>
<td><img width=129 src="TreeK-D=3-R=2-L=4-y-1.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=5-y-1.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=6-y-1.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=7-y-1.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=8-y-1.png"></img></td>
</tr>
<tr>
<th>project = Z</th>
<td><img width=129 src="TreeK-D=3-R=2-L=4-z-1.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=5-z-1.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=6-z-1.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=7-z-1.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=2-L=8-z-1.png"></img></td>
</tr>
</table>


<h3>TreeK-D3-R4-L?</h3>

<?php tests("Mesh","test_TreeK","test_TreeK-D3-R4-L4"); ?>
<?php tests("Mesh","test_TreeK","test_TreeK-D3-R4-L6"); ?>
<?php tests("Mesh","test_TreeK","test_TreeK-D3-R4-L8"); ?>

<table>
<tr>
<th></th>
<th colspan=3>coalesce = false</th>
<th colspan=3>coalesce = true</th>
</tr>
<tr>
<th></th>
<th>levels = 4</th>
<th>levels = 6</th>
<th>levels = 8</th>
<th>levels = 4</th>
<th>levels = 6</th>
<th>levels = 8</th>
</tr>
<tr>
<th>project = X</th>
<td><img width=129 src="TreeK-D=3-R=4-L=4-x-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=4-L=6-x-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=4-L=8-x-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=4-L=4-x-1.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=4-L=6-x-1.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=4-L=8-x-1.png"></img></td>
</tr>
<tr>
<th>project =  Y</th>
<td><img width=129 src="TreeK-D=3-R=4-L=4-y-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=4-L=6-y-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=4-L=8-y-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=4-L=4-y-1.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=4-L=6-y-1.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=4-L=8-y-1.png"></img></td>
</tr>
<tr>
<th>project = Z</th>
<td><img width=129 src="TreeK-D=3-R=4-L=4-z-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=4-L=6-z-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=4-L=8-z-0.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=4-L=4-z-1.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=4-L=6-z-1.png"></img></td>
<td><img width=129 src="TreeK-D=3-R=4-L=8-z-1.png"></img></td>
</br/>
</body>
</html>
