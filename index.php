<html>
  <head>
    <title>Cello Test Page</title>
    <link href="cello.css" rel="stylesheet" type="text/css">
  </head>
<body>
  <h1>Cello Test Page</h1>

<?php
function pass($testrun) {
     echo '<span style="color:blue">';
     system ("awk '/pass/{print \"passed: \",$5,\"</br/>\"}' < test/test_$testrun.out");
     echo '</span>';
}

function fail($testrun) {
     echo '<span style="color:red">';
     system ("awk '/FAIL/{print \"FAILED: \",$3,\"</br/>\"}' < test/test_$testrun.out");
     echo '</span>';
}
function tests($testrun) {
  pass($testrun);
  fail($testrun);
}
?>


<h3>Data</h3>
<?php tests("Data"); ?>
<?php tests("FieldBlock"); ?>
<?php tests("FieldDescr"); ?>
<h2>Disk</h2>
<?php tests("FileHdf5"); ?>
<?php tests("Ifrit"); ?>
<h2>Distribute</h2>
<h2>Error</h2>
<?php tests("Error"); ?>
<h2>Field</h2>
<h2>Global</h2>
<h2>Memory</h2>
<?php tests("Memory"); ?>
<h2>Mesh</h2>
<h2>Method</h2>
<h2>Monitor</h2>
<?php tests("Monitor"); ?>
<img src="test1.png"></img>
<img src="test2.png"></img>
<img src="test3.png"></img>
<img src="test4.png"></img>

<h2>Parallel</h2>
<?php tests("parallel"); ?>
<h2>Parameters</h2>
<?php tests("Parameters"); ?>
<h2>Particles</h2>
<h2>Performance</h2>
<?php tests("Performance"); ?>
<h2>Portal</h2>
<h2>Schedule</h2>
<?php tests("Schedule"); ?>

<h2>Simulation</h2>
<?php tests("Simulation"); ?>
<h2>Task</h2>
<h2>Test</h2>
<h2>User</h2>


<h4>ppm_image</h4>

<img width=320 src="slice-ppm-image-000000.png"></img>
<img width=320  src="slice-ppm-image-000010.png"></img>
<h4>ppm_implosion</h4>
<img width=320 src="slice-ppm-implosion-000000.png"></img>
<img width=320 src="slice-ppm-implosion-000010.png"></img>
<table>

<h4>ppm_implosion3</h4>

<img width=200 src="project-ppm-implosion3-000000-x.png"></img>
<img width=200 src="project-ppm-implosion3-000000-y.png"></img>
<img width=200 src="project-ppm-implosion3-000000-z.png"></img>
<img width=200 src="project-ppm-implosion3-000010-x.png"></img>
<img width=200 src="project-ppm-implosion3-000010-y.png"></img>
<img width=200 src="project-ppm-implosion3-000010-z.png"></img>

<h3>ppml_blast</h3>

<img width=200 src="project-ppml-blast-000000-x.png"></img>
<img width=200 src="project-ppml-blast-000000-y.png"></img>
<img width=200 src="project-ppml-blast-000000-z.png"></img>
<img width=200 src="project-ppml-blast-000010-x.png"></img>
<img width=200 src="project-ppml-blast-000010-y.png"></img>
<img width=200 src="project-ppml-blast-000010-z.png"></img>

<h3>ppml_implosion</h3>

<img width=200 src="project-ppml-implosion3-000000-x.png"></img>
<img width=200 src="project-ppml-implosion3-000000-y.png"></img>
<img width=200 src="project-ppml-implosion3-000000-z.png"></img>
<img width=200 src="project-ppml-implosion3-000010-x.png"></img>
<img width=200 src="project-ppml-implosion3-000010-y.png"></img>
<img width=200 src="project-ppml-implosion3-000010-z.png"></img>

<h3>MethodEnzoPpm</h3>
<tr>
<th>density</th>
<th>velocity_x</th>
<th>velocity_y</th>
<th>total_energy</th></tr>
<tr>
<td><img width=200 src="ppm-density-00000.png"></img></td>
<td><img width=200 src="ppm-velocity_x-00000.png"></img></td>
<td><img width=200 src="ppm-velocity_y-00000.png"></img></td>
<td><img width=200 src="ppm-total_energy-00000.png"></img></td>
</tr>
<tr>
<td><img width=200 src="ppm-density-00100.png"></img></td>
<td><img width=200 src="ppm-velocity_x-00100.png"></img></td>
<td><img width=200 src="ppm-velocity_y-00100.png"></img></td>
<td><img width=200 src="ppm-total_energy-00100.png"></img></td>
</tr>
</table>

<h3>TreeK-D2-R2-L?</h3>
<img width=257 src="TreeK-D=2-R=2-L=6-0.png"></img>
<img width=257 src="TreeK-D=2-R=2-L=7-0.png"></img>
<img width=257 src="TreeK-D=2-R=2-L=8-0.png"></img>
<img width=257 src="TreeK-D=2-R=2-L=9-0.png"></img>
<img width=257 src="TreeK-D=2-R=2-L=10-0.png"></img>
</br/>
<img width=257 src="TreeK-D=2-R=2-L=6-1.png"></img>
<img width=257 src="TreeK-D=2-R=2-L=7-1.png"></img>
<img width=257 src="TreeK-D=2-R=2-L=8-1.png"></img>
<img width=257 src="TreeK-D=2-R=2-L=9-1.png"></img>
<img width=257 src="TreeK-D=2-R=2-L=10-1.png"></img>

<h3>TreeK-D2-R4-L?</h3>
<img width=257 src="TreeK-D=2-R=4-L=6-0.png"></img>
<img width=257 src="TreeK-D=2-R=4-L=8-0.png"></img>
<img width=257 src="TreeK-D=2-R=4-L=10-0.png"></img>
</br/>
<img width=257 src="TreeK-D=2-R=4-L=6-1.png"></img>
<img width=257 src="TreeK-D=2-R=4-L=8-1.png"></img>
<img width=257 src="TreeK-D=2-R=4-L=10-1.png"></img>

<h3>TreeK-D3-R2-L?</h3>

<img width=129 src="TreeK-D=3-R=2-L=4-x-0.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=5-x-0.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=6-x-0.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=7-x-0.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=8-x-0.png"></img>
</br/>
<img width=129 src="TreeK-D=3-R=2-L=4-y-0.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=5-y-0.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=6-y-0.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=7-y-0.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=8-y-0.png"></img>
</br/>
<img width=129 src="TreeK-D=3-R=2-L=4-z-0.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=5-z-0.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=6-z-0.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=7-z-0.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=8-z-0.png"></img>
</br/>

<img width=129 src="TreeK-D=3-R=2-L=4-x-1.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=5-x-1.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=6-x-1.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=7-x-1.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=8-x-1.png"></img>
</br/>
<img width=129 src="TreeK-D=3-R=2-L=4-y-1.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=5-y-1.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=6-y-1.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=7-y-1.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=8-y-1.png"></img>
</br/>
<img width=129 src="TreeK-D=3-R=2-L=4-z-1.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=5-z-1.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=6-z-1.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=7-z-1.png"></img>
<img width=129 src="TreeK-D=3-R=2-L=8-z-1.png"></img>
</br/>

<h3>TreeK-D3-R4-L?</h3>
<img width=129 src="TreeK-D=3-R=4-L=4-x-0.png"></img>
<img width=129 src="TreeK-D=3-R=4-L=6-x-0.png"></img>
<img width=129 src="TreeK-D=3-R=4-L=8-x-0.png"></img>
<img width=129 src="TreeK-D=3-R=4-L=4-x-1.png"></img>
<img width=129 src="TreeK-D=3-R=4-L=6-x-1.png"></img>
<img width=129 src="TreeK-D=3-R=4-L=8-x-1.png"></img>
</br/>
<img width=129 src="TreeK-D=3-R=4-L=4-y-0.png"></img>
<img width=129 src="TreeK-D=3-R=4-L=6-y-0.png"></img>
<img width=129 src="TreeK-D=3-R=4-L=8-y-0.png"></img>
<img width=129 src="TreeK-D=3-R=4-L=4-y-1.png"></img>
<img width=129 src="TreeK-D=3-R=4-L=6-y-1.png"></img>
<img width=129 src="TreeK-D=3-R=4-L=8-y-1.png"></img>
</br/>
<img width=129 src="TreeK-D=3-R=4-L=4-z-0.png"></img>
<img width=129 src="TreeK-D=3-R=4-L=6-z-0.png"></img>
<img width=129 src="TreeK-D=3-R=4-L=8-z-0.png"></img>
<img width=129 src="TreeK-D=3-R=4-L=4-z-1.png"></img>
<img width=129 src="TreeK-D=3-R=4-L=6-z-1.png"></img>
<img width=129 src="TreeK-D=3-R=4-L=8-z-1.png"></img>
</br/>
</body>
</html>
