<html>
  <head>
    <title>Cello Unit Test Page</title>
    <link href="cello.css" rel="stylesheet" type="text/css">
  </head>
<body>
  <h1>Cello Unit Test Page</h1>

<?php
function tests($component,$testrun) {
  
   $file = "src/$component/test_$testrun.cpp";
   $output = "test/test_$testrun.unit";

   echo "</br/>\n";
   echo "<table><tr>\n";
   echo "<td ><a href=\"$file\">$file</a> &gt; <a href=\"$output\">$output</a></td></tr>\n";
   echo "</table>\n";
   echo "</br/>\n";
   system("awk 'BEGIN{c=0}; /UNIT TEST END/ {c=1}; END{ if (c==0) print \"<strong class=fail>TEST PROGRAM DID NOT COMPLETE\!</strong>\"; if (c!=0) print \"<strong class=pass>Test program completed normally</strong>\"}' < $output");
   echo "</br/>\n";
   echo "</br/>\n";
   test($testrun,"FAIL",$component,$file);
   echo "</br/>\n";
   test($testrun,"pass",$component,$file);
 };

function test($testrun,$type,$component) {

  if ($type == "pass") {
    $cols = "$4,$6";
    $itemtext  = "";
    $rowtext = "</tr><tr>";
  } else {
    $cols = "\"$file\",$4,$6";
    $itemtext  = "</tr><tr>";
    $rowtext = "";
  }

  $output = "test/test_$testrun.unit";
  echo "${type}ed tests: </br/>";
  echo "<table><tr>";
  system ("awk 'BEGIN {c=0}; / $type /{split($3,a,\"\/\"); print \"<td class=$type> \",$cols , \" </td>$itemtext\"; c=c+1}; {if (c==5) {c=0; print \"$rowtext\"}}' < $output");
  echo "</tr></table>";
     
};

?>

<p>
Test data shown on this page is automatically generated whenever <code>Cello</code> is compiled.  FAIL status usually indicates that the corresponding function has not been implemented yet.  Empty tables indicate that the unit test files have been deleted, probably with "<code>scons -c</code>" (<code>scons</code> version of "<code>make clean</code>").
</p>

<h3>Data</h3>

<?php tests("Data","Data"); ?>

<h2>Disk</h2>

<?php tests("Disk","FileHdf5"); ?>
<?php tests("Disk","Ifrit"); ?>

<h2>Distribute</h2>




<h2>Error</h2>

<?php tests("Error","Error"); ?>

<h2>Field</h2>

<?php tests("Field","FieldBlock"); ?>
<!-- <?php tests("FieldBlockFaces"); ?> -->
<?php tests("Field","FieldDescr"); ?>

<h2>Global</h2>


<h2>Memory</h2>

<?php tests("Memory","Memory"); ?>

<h2>Mesh</h2>


<h2>Method</h2>


<h2>Monitor</h2>

<?php tests("Monitor","Monitor"); ?>
<img src="monitor_image_1.png"></img>
<img src="monitor_image_2.png"></img>
<img src="monitor_image_3.png"></img>
<img src="monitor_image_4.png"></img>


<h2>Parallel</h2>

<?php tests("Parallel","parallel"); ?>

<h2>Parameters</h2>

<?php tests("Parameters","Parameters"); ?>

<h2>Particles</h2>


<h2>Performance</h2>

<?php tests("Performance","Performance"); ?>

<h2>Portal</h2>


<h2>Schedule</h2>

<?php tests("Schedule","Schedule"); ?>


<h2>Simulation</h2>

<?php tests("Simulation","Simulation"); ?>

<h2>Task</h2>


<h2>Test</h2>

<h2>User</h2>


<h4>test_ppm ppm_image</h4>

<p>TODO: Add separate program and test output files</p>

<?php tests("Enzo","ppm_image"); ?>

<table>
<tr><th>cycle = 0</th><th>cycle=10</th></tr>
<tr>
<td><img width=320 src="slice-ppm-image-000000.png"></img></td>
<td><img width=320  src="slice-ppm-image-000010.png"></img></td>
</tr>
</table>

<h4>test_ppm ppm_implosion</h4>

<table>
<tr><th>cycle = 0</th><th>cycle=10</th></tr>
<tr>
<td><img width=320 src="slice-ppm-implosion-000000.png"></img></td>
<td><img width=320 src="slice-ppm-implosion-000010.png"></img></td>
</tr>
</table>

<h4>test_ppm ppm_implosion3</h4>

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

<h3>MethodEnzoPpm</h3>

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
