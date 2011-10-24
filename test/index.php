<html>
<head>
<title>Enzo-P / Cello Test Results</title>
<link href="cello.css" rel="stylesheet" type="text/css">
<?php
if (file_exists("COMPILING")) {
  echo "<meta http-equiv=\"refresh\" content=60>";
}
?>
   </head>
   <body>
   <h1>Enzo-P / Cello Test Results</h1>

   <?php

   //----------------------------------------------------------------------

   function test_group($testgroup) {
   printf ("<hr><a name=\"$testgroup\"><h2>$testgroup Tests</h2>");
 }

   //----------------------------------------------------------------------

   function test_output ($output_file) {
     if (file_exists($output_file)) {
       $output_html = "<a href=\"$output_file\">$output.unit</a>";
       system("cat $output_file | awk 'BEGIN{c=0}; /UNIT TEST END/ {c=1}; END{ if (c==0) print \"<td class=fail><a href='$output_file'>incomplete</a></td>\"; if (c!=0) print \"<td class=pass><a href='$output_file'>complete</a></td>\"}'");
     } else {
       echo "<td></td>";
     }
   }

   //----------------------------------------------------------------------

function test_date ($output_file) {
  if (file_exists($output_file)) {
    $output_html = date ("Y-m-d", filemtime($output_file));
    echo "<td class=pass>$output_html</td>";
  } else {
    echo "<td class=pass></td>";
  }
}

   //----------------------------------------------------------------------

function test_time ($output_file) {
  if (file_exists($output_file)) {
    $output_html = date ("H:i:s", filemtime($output_file));
    echo "<td class=pass>$output_html</td>";
  } else {
    echo "<td class=pass></td>";
  }
}  

   //----------------------------------------------------------------------

function test_passed ($output_file) {
  if (file_exists($output_file)) {
    echo "<td class=pass>";
    system("cat $output_file | grep pass | grep '0/' | wc -l");
    echo "</td>";
  } else {
    echo "<td class=pass></td>";
  }
}

   //----------------------------------------------------------------------

function test_unfinished ($output_file) {
  if (file_exists($output_file)) {
    system("grep 0/ $output_file | awk 'BEGIN{c=0}; /incomplete/ {c=c+1}; END {if (c!=0) print\"<td class=unfinished>\",c,\"</td>\"; if (c==0) print \"<td class=pass>0</td>\"}'");
  } else {
    echo "<td class=pass></td>";
  }
}

   //----------------------------------------------------------------------

function test_failed ($output_file) {
  if (file_exists($output_file)) {
    system("grep 0/ $output_file | awk 'BEGIN{c=0}; /FAIL/ {c=c+1}; END {if (c!=0) print\"<td class=fail>\",c,\"</td>\"; if (c==0) print \"<td class=pass>0</td>\"}'");
  } else {
    echo "<td class=pass></td>";
  }
}

   //----------------------------------------------------------------------

$types = array("mpi","charm");
  $num_types = sizeof($types);

  function tests($component,$testrun,$output,$test_name) {

  global $types;
  global $num_types;

  $input_file  = "../input/$output.in";
  $input_link  = "../cello-src/input/$output.in";
  $input_html = "<a href=\"$input_link\">$output.in</a>";

  $source_file = "../src/$component/$testrun.cpp";
  $source_link = "../cello-src/src/$component/$testrun.cpp";
  $source_html = "<a href=\"$source_link\">$testrun.cpp</a>";

  if (! file_exists($source_file)) {
    echo "<p><strong class=fail>$source_file does not exist</strong></p>";
  } else {
    echo "<table>\n";

    echo "<tr><th colspan=7>$source_html \n";
    if ($test_name != "") {
      echo "($test_name)";
    }
    if (file_exists($input_file)) {
      echo " $input_html";
    }
    echo "</th></tr>\n";
    echo "<tr>\n";
    echo "   <th></th>";
    echo "   <th>Output</th>";
    echo "   <th>Date</th>";
    echo "   <th>Time</th>";
    echo "   <th>Failed</th>";
    echo "   <th>Unfinished</th>";
    echo "   <th>Passed</th>";
    echo "</tr>\n";

    //--------------------------------------------------

    for ($i = 0; $i<$num_types; ++$i) {

      echo "<tr>\n";

      echo "<th> $types[$i] </th>";

      $output_file = "$types[$i]/$output.unit";

      test_output     ($output_file);
      test_date       ($output_file);
      test_time       ($output_file);
      test_failed     ($output_file);
      test_unfinished ($output_file);
      test_passed     ($output_file);

      echo "</tr>\n";
    }

    echo "</tr></table></br/>\n";

    echo "<table><tr>";

    for ($i = 0; $i<$num_types; ++$i) {
      $type = $types[$i];
      $output_file = "../test/$type/$output.unit";
      if (file_exists($output_file)) {
	test($type,$testrun,"FAIL");
      }
    }

    echo "</tr></table></br/>";
  }
};

   //----------------------------------------------------------------------

function test($type,$testrun,$type) {
  $ltype = strtolower($type);

  $cols = "$4,$6,$7,$8,$9,$10";
  $rowtext = "</tr><tr>";

  $output = "../test/$type/$testrun.unit";
  $count = exec("cat $output | grep $type | grep '0/' | wc -l");
  if ($count == 0) {
    #     echo "<strong >no ${ltype}ed tests</strong></br/>";
  } else {
    echo "<th class=$type><strong>$type ${ltype}ed</strong></th> ";
    system ("grep '0/' $output | sort | uniq | awk 'BEGIN {c=1}; / $type /{split($3,a,\"\/\"); print \"<td class=$type> \",$cols , \" </td>\"; c=c+1}; {if (c==5) {c=0; print \"$rowtext\"}}'");
    echo "</tr><tr></tr>";
  }
     
};

?>

This page contains the current status of Enzo-P / Cello unit tests, as run
on the main development platform.

  <hr>

  <h2>Test Summary</h2>

  <?php

   //----------------------------------------------------------------------

  function summary_missing_executable ($test_output, $executables)
{
  global $types;
  global $num_types;

    for ($i = 0; $i<$num_types; ++$i) {

    $type = $types[$i];

    $count_missing = 0;

    for ($test = 0; $test<sizeof($executables); ++$test) {
      $exe = $executables[$test];
      $bin = "../bin/$type/$exe";
      if (! file_exists($bin)) {
        ++ $count_missing ;
      }
    }

    if ($count_missing == 0) {
      printf ("<td></td>");
    } else {
      printf ("<td class=fail>$count_missing</td>");
    }
  }
  printf ("<th></th>");
}

   //----------------------------------------------------------------------

function summary_missing_output ($test_output, $executables)
{
  global $types;
  global $num_types;

  for ($i = 0; $i<$num_types; ++$i) {

    $count_missing = 0;
    for ($test = 0; $test<sizeof($test_output); ++$test) {
      $output = "../test/$types[$i]/test_$test_output[$test].unit";
      if (! file_exists($output)) {
	++ $count_missing;
      }
    }
    if ($count_missing == 0) {
      printf ("<td></td>");
    } else {
      printf ("<td class=fail>$count_missing</td>");
    }
  }
  printf ("<th></th>");
}

   //----------------------------------------------------------------------

function summary_crashed_runs ( $test_output, $executables)
{
  global $types;
  global $num_types;

  for ($i = 0; $i<$num_types; ++$i) {

    $output_files = "";
    $num_output_files = 0;
    for ($test = 0; $test<sizeof($test_output); ++$test) {
      $output = "../test/$types[$i]/test_$test_output[$test].unit";
      $output_files = "$output_files $output";
      ++$num_output_files;
    }
    system("cat $output_files | awk 'BEGIN{b=0; e=0;}; /UNIT TEST BEGIN/ {b=b+1};/UNIT TEST END/ {e=e+1};END{if (b==e) {print \"<td></td>\"} else {print \"<td class=fail>\"b - e\"</td>\";}}'");
  }
  printf ("<th></th>");
}

   //----------------------------------------------------------------------

function summary_failed_tests ($test_output, $executables)
{
  global $types;
  global $num_types;

  for ($i = 0; $i<$num_types; ++$i) {

    $output_files = "";
    for ($test = 0; $test<sizeof($test_output); ++$test) {
      $output = "../test/$types[$i]/test_$test_output[$test].unit";
      $output_files = "$output_files $output";
    }

    system("grep '0/' $output_files | awk 'BEGIN {c=0}; /FAIL/{c=c+1}; END{if (c==0) {print \"<td></td>\"} else {print \"<td class=fail>\",c,\"</td>\";}} '");
  }
  printf ("<th></th>");
}

   //----------------------------------------------------------------------

function summary_unfinished_tests ($test_output, $executables)
{
  global $types;
  global $num_types;

  for ($i = 0; $i<$num_types; ++$i) {

    $output_files = "";
    for ($test = 0; $test<sizeof($test_output); ++$test) {
      $output = "../test/$types[$i]/test_$test_output[$test].unit";
      $output_files = "$output_files $output";
    }

    system("grep '0/' $output_files | awk 'BEGIN {c=0}; /incomplete/{c=c+1}; END{if (c==0) {print \"<td></td>\"} else {print \"<td class=unfinished>\",c,\"</td>\";}} '");
  }
  printf ("<th></th>");
}

   //----------------------------------------------------------------------

function summary_passed_tests ($test_output, $executables)
{

  global $types;
  global $num_types;

  for ($i = 0; $i<$num_types; ++$i) {

    $output_files = "";
    for ($test = 0; $test<sizeof($test_output); ++$test) {
      $output = $test_output[$test];
      $output_files = "$output_files ../test/$types[$i]/test_$output.unit";
    }
    system("grep '0/' $output_files | awk 'BEGIN {c=0}; /pass/{c=c+1}; END{if (c==0) {print \"<td></td>\"} else {print \"<td class=pass>\",c,\"</td>\";}} '");

  }
  printf ("<th></th>");
}

   //----------------------------------------------------------------------

  function test_summary($component,$test_output,$executables)
{
  printf ("<tr><th><a href=\"#$component\">$component</a></th>\n");

  global $types;

  // Missing executable

  summary_missing_executable ($test_output, $executables);
  summary_crashed_runs       ($test_output, $executables);
  summary_failed_tests       ($test_output, $executables);
  summary_unfinished_tests   ($test_output, $executables);
  summary_passed_tests       ($test_output, $executables);

  printf ("</tr>\n");
}

function test_table ($file_root,$size_array, $types)
{
  echo "<table>";
  echo "<tr>";
    echo "<th>$file_root</th>  <th>movie</th>";
    for ($j = 0; $j < sizeof($size_array); ++$j) {
      $size = $size_array[$j];
      printf ("<th>$size</th>\n");
    }
    echo "</tr>";
    for ($i = 0; $i < sizeof($types); ++$i) {
      echo "<tr>";
      $type = $types[$i];
     	printf ("<th>$type</th>\n"); 
	// Show movie file if available
	echo "<td>";
         $swf_file = "$type/$file_root.swf"; 
	 $size_last = $size_array[sizeof($size_array)-1]; 
	 $png_file_last = "$type/$file_root-$size_last.png"; 
	 swf_movie($swf_file, 
     	   	  $png_file_last, 
     	   	  160); 
	echo "</td>";
	// Show available image frames
	for ($j = 0; $j < sizeof($size_array); ++$j) {
	  $size = $size_array[$j];
          $png_file = "$type/$file_root-$size.png"; 

     	   printf ("<td><img width=160 src=$png_file></img></td>\n");  
     	 }  
     	 echo "</tr>";  
    }
    echo "</table></br>";
}
  function swf_movie ($filename, $last_image, $image_size)
  {
    global $types;
    global $num_types;
if (file_exists($last_image)) {
  printf ("<OBJECT classid=\"clsid:D27CDB6E-AE6D-11cf-96B8-444553540000\"\n");
  printf ("codebase=\"http://download.macromedia.com/pub/shockwave/cabs/flash/swflash.cab#version=6,0,0,0\"\n");
  printf ("WIDTH=\"$image_size\" HEIGHT=\"$image_size\"\n");
  printf ("        id=\"implosion\" ALIGN=\"\">\n");
printf ("     <PARAM NAME=$filename\n");
printf ("            VALUE=\"$filename\">\n");
printf ("     <PARAM NAME=quality VALUE=high>\n");
printf ("     <PARAM NAME=bgcolor VALUE=#333399>\n");
printf ("     <EMBED src=\"$filename\"\n");
printf ("          quality=high\n");
printf ("          bgcolor=#333399\n");
printf ("          WIDTH=\"$image_size\" HEIGHT=\"$image_size\"\n");
printf ("          NAME=\"$filename\" ALIGN=\"\"\n");
printf ("          TYPE=\"application/x-shockwave-flash\"\n");
printf ("          PLUGINSPAGE=\"http://www.macromedia.com/go/getflashplayer\">\n");
printf ("       </EMBED>\n");
printf ("      </OBJECT> \n");
}  
  }

printf ("<table>\n");
printf ("<tr>\n");

  if (file_exists("COMPILING"))  {
    printf ( "<th rowspan=2 class=compiling>");
    printf ("<strong> COMPILING </strong>\n");
    printf ("</th>");
  } else { 
    printf ("<th rowspan=2></th>\n"); 
  }
  
printf ( "<th colspan=$num_types class=fail>Missing Executable</th>");
printf ("<th></th>");
// printf ( "<th colspan=$num_types class=fail>Missing Output</th>");
// printf ("<th></th>");
printf ( "<th colspan=$num_types class=fail>Crashed Runs</th>");
printf ("<th></th>");
printf ( "<th colspan=$num_types class=fail>Failed Tests</th>");
printf ("<th></th>");
printf ( "<th colspan=$num_types class=unfinished>Unfinished Tests</th>");
printf ("<th></th>");
printf ( "<th colspan=$num_types class=pass>Passed Tests</th>");
printf ("<th></th>");
printf ( "</tr><tr>\n");

for ($k = 0; $k < 5; $k ++) {
  for ($i = 0; $i < $num_types; ++$i) {
     $type_active = "";
     if (file_exists("COMPILING"))  {
        $type_active = file_get_contents("COMPILING");
     }
     if ($type_active == $types[$i]) {
        printf ("<th class=compiling>");
     } else {
        printf ("<th> ");
     }
     printf (" <a href=$types[$i]/out.scons>$types[$i]</a> </th>");
  }
  printf ("<th> </th>");
}


test_summary("Enzo-PPM",
	     array("enzo-p_1",
		   "enzo-p_2"),
	     array("enzo-p",  "enzo-p"));

test_summary("Enzo-BC-2D", 
	     array("boundary_reflecting-2d",
		   "boundary_periodic-2d",
		   "boundary_inflow-2d",
		   "boundary_outflow-2d"),
	     array("enzo-p", "enzo-p", "enzo-p", "enzo-p"));

test_summary("Enzo-BC-3D",
	     array("boundary_reflecting-3d",
		   "boundary_periodic-3d",
		   "boundary_inflow-3d",
		   "boundary_outflow-3d"),
	     array("enzo-p", "enzo-p", "enzo-p", "enzo-p"));

// Print row divider

printf ("<tr><th></th>");
for ($k = 0; $k < 5; $k ++) {
  for ($i = 0; $i < $num_types; ++$i) {
    printf ("<th> </th>");
  }
  printf ("<th></th>");
}
printf ("</tr>\n");


test_summary("Disk",array(     "FileHdf5",     "FileIfrit"),
	     array("test_FileHdf5","test_FileIfrit")); 
test_summary("Error",array(    "Error"),
	     array("test_Error")); 
test_summary("Field",array(    "FieldBlock",     "FieldDescr",     "FieldFace",     "ItField"),
	     array("test_FieldBlock","test_FieldDescr","test_FieldFace","test_ItField")); 
test_summary("Io",array("ItReduce"),
	     array("test_ItReduce")); 
test_summary("Memory",array("Memory"),
	     array("test_Memory")); 
test_summary("Mesh",array("Hierarchy","Patch","Block"),
	     array("test_Hierarchy","test_Patch","test_Block")); 
test_summary("Monitor",array("Monitor"),
	     array("test_Monitor")); 
test_summary("Parallel",array("GroupProcess","Layout"),
	     array("test_GroupProcess","test_Layout")); 
test_summary("Parameters",array("Parameters"),
	     array("test_Parameters")); 
test_summary("Performance",array("Papi", "Performance"),
	     array("test_Papi","test_Performance")); 


printf ("</tr></table></br>\n");

//======================================================================

test_group("Enzo-PPM");

echo "<h3>Implosion2D (serial) </h3>";

?>

Enzo-PPM tests serve to test basic PPM functionality in Enzo-P.  A
small implosion problem is run for 400 cycles, first with
one block (1,1) then eight blocks (2,4).

</p>

Currently, "serial" results are incorrect for multiple blocks, which
is to be expected.  There are errors in parallel CHARM++ and MPI with
eight blocks because the final time after 400 cycles does not exactly
match the time for the serial runs.  The results look qualitatively
correct however, even at time 2.5 for 400<sup>2</sup>(over 13000
cycles).
</p>

<?php

tests("Enzo","enzo-p","test_enzo-p_1","PPM 1 block");

test_table ("enzo-p_1-d",
	    array("000000","000100","000200"), $types);

//----------------------------------------------------------------------

echo "<h3>Implosion2D (parallel) </h3>";

tests("Enzo","enzo-p","test_enzo-p_2","PPM 8 blocks");

test_table ("enzo-p_2-d",
	    array("000000","000100","000200"), $types);


//======================================================================

  test_group("Enzo-BC-2D");

echo "<h3>2D Reflecting</h3>";

tests("Enzo","enzo-p","test_boundary_reflecting-2d","Reflecting 2D");

test_table ("boundary_reflecting-2d",
	    array("0000","0100","0200","0300","0400","0500"), $types);

//----------------------------------------------------------------------

echo "<h3>2D Periodic</h3>";

tests("Enzo","enzo-p","test_boundary_periodic-2d","Periodic 2D");

test_table ("boundary_periodic-2d",
	    array("0000","0100","0200","0300","0400","0500"), $types);

//----------------------------------------------------------------------

echo "<h3>2D Inflow</h3>";

tests("Enzo","enzo-p","test_boundary_inflow-2d","Inflow 2D");

test_table ("boundary_inflow-2d",
	    array("0000","0100","0200","0300","0400","0500"), $types);

//----------------------------------------------------------------------

echo "<h3>2D Outflow</h3>";
  
tests("Enzo","enzo-p","test_boundary_outflow-2d","Outflow 2D");

test_table ("boundary_outflow-2d",
	    array("0000","0100","0200","0300","0400","0500"), $types);

//======================================================================

test_group("Enzo-BC-3D");


echo "<h3>3D Reflecting</h3>";

tests("Enzo","enzo-p","test_boundary_reflecting-3d","Reflecting 3D");

test_table ("boundary_reflecting-3d",
	    array("0000","0020","0040","0060","0080"), $types);

//----------------------------------------------------------------------

echo "<h3>3D Periodic</h3>";

tests("Enzo", "enzo-p","test_boundary_periodic-3d","Periodic 3D");

test_table ("boundary_periodic-3d",
	    array("0000","0020","0040","0060","0080"), $types);

//----------------------------------------------------------------------

echo "<h3>3D Inflow</h3>";

tests("Enzo","enzo-p","test_boundary_inflow-3d","Inflow 3D");

test_table ("boundary_inflow-3d",
	    array("0000","0020","0040","0060","0080"), $types);

//----------------------------------------------------------------------

echo "<h3>3D Outflow</h3>";

tests("Enzo","enzo-p","test_boundary_outflow-3d","Outflow 3D");

test_table ("boundary_outflow-3d",
	    array("0000","0020","0040","0060","0080"), $types);

   //======================================================================

test_group("Disk");

tests("Disk","test_FileHdf5", "test_FileHdf5","");
tests("Disk","test_FileIfrit","test_FileIfrit","");

   //----------------------------------------------------------------------

test_group("Error");

tests("Error","test_Error","test_Error","");


   //----------------------------------------------------------------------

test_group("Field");

tests("Field","test_FieldDescr","test_FieldDescr","");
tests("Field","test_FieldBlock","test_FieldBlock","");
tests("Field","test_FieldFace","test_FieldFace","");
tests("Field","test_ItField","test_ItField","");

   //----------------------------------------------------------------------

test_group("Io");

tests("Io","test_ItReduce", "test_ItReduce","");

   //----------------------------------------------------------------------

test_group("Memory");

tests("Memory","test_Memory","test_Memory","");


   //----------------------------------------------------------------------

test_group("Mesh");

tests("Mesh","test_Hierarchy","test_Hierarchy",""); 
tests("Mesh","test_Patch","test_Patch",""); 
tests("Mesh","test_Block","test_Block",""); 

   //----------------------------------------------------------------------

test_group("Monitor");

tests("Monitor","test_Monitor","test_Monitor","");

// printf ("<img src=\"monitor_image_1.png\"></img>\n");
// printf ("<img src=\"monitor_image_2.png\"></img>\n");
// printf ("<img src=\"monitor_image_3.png\"></img>\n");
// printf ("<img src=\"monitor_image_4.png\"></img>\n");

   //----------------------------------------------------------------------

test_group("Parallel");

tests("Parallel","test_GroupProcess","test_GroupProcess","");
tests("Parallel","test_Layout","test_Layout","");

   //----------------------------------------------------------------------

test_group("Parameters");

tests("Parameters","test_Parameters","test_Parameters","");

   //----------------------------------------------------------------------

test_group("Performance");

tests("Performance","test_Performance","test_Performance","");
tests("Performance","test_Papi",       "test_Papi","");
?>



  <h3>TreeK-D2-R2-L?</h3>

  <?php
  tests("Mesh","test_TreeK","test_TreeK-D2-R2-L6", "2D L=6 r=2");
  tests("Mesh","test_TreeK","test_TreeK-D2-R2-L7", "2D L=7 r=2");
  tests("Mesh","test_TreeK","test_TreeK-D2-R2-L8", "2D L=8 r=2");
  ?>

<table>
<tr>
<th>coalesce</th>
<th>levels = 6</th>
  <th>levels = 7</th>
  <th>levels = 8</th>
  </tr>
  <tr>
  <th>false</th>
  <td><img width=257 src="serial/TreeK-D=2-R=2-L=6-0.png"></img></td>
  <td><img width=257 src="serial/TreeK-D=2-R=2-L=7-0.png"></img></td>
  <td><img width=257 src="serial/TreeK-D=2-R=2-L=8-0.png"></img></td>
  </tr>
  <tr>
  <th>true</th>
  <td><img width=257 src="serial/TreeK-D=2-R=2-L=6-1.png"></img></td>
  <td><img width=257 src="serial/TreeK-D=2-R=2-L=7-1.png"></img></td>
  <td><img width=257 src="serial/TreeK-D=2-R=2-L=8-1.png"></img></td>
  </tr>
  </table></br>

  <h3>TreeK-D2-R4-L?</h3>

  <?php
  tests("Mesh","test_TreeK","test_TreeK-D2-R4-L6", "2D L=6 r=4");
  tests("Mesh","test_TreeK","test_TreeK-D2-R4-L8", "2D L=8 r=4");
  ?>

<table>
<tr>
<th>coalesce</th>
<th>levels = 6</th>
  <th>levels = 8</th>
  <!-- <th>levels = 10</th> -->
  </tr>
  <tr>
  <th>false</th>
  <td><img width=257 src="serial/TreeK-D=2-R=4-L=6-0.png"></img></td>
  <td><img width=257 src="serial/TreeK-D=2-R=4-L=8-0.png"></img></td>
  <!-- <td><img width=257 src="serial/TreeK-D=2-R=4-L=10-0.png"></img></td> -->
  </tr>
  <tr>
  <th>true</th>
  <td><img width=257 src="serial/TreeK-D=2-R=4-L=6-1.png"></img></td>
  <td><img width=257 src="serial/TreeK-D=2-R=4-L=8-1.png"></img></td>
  <!-- <td><img width=257 src="serial/TreeK-D=2-R=4-L=10-1.png"></img></td> -->
  </tr>
  </table></br>

  <h3>TreeK-D3-R2-L?</h3>

  <?php
  tests("Mesh","test_TreeK","test_TreeK-D3-R2-L4", "3D L=4 r=2");
  tests("Mesh","test_TreeK","test_TreeK-D3-R2-L5", "3D L=5 r=2");
  tests("Mesh","test_TreeK","test_TreeK-D3-R2-L6", "3D L=6 r=2");
   ?>

<table>
<tr>
<th>coalesce = false</th>
  <th>levels = 4</th>
  <th>levels = 5</th>
  <th>levels = 6</th>
  <!-- <th>levels = 7</th> -->
  <!-- <th>levels = 8</th> -->
  </tr>
  <tr>
  <th>project = X</th>
  <td><img width=129 src="serial/TreeK-D=3-R=2-L=4-x-0.png"></img></td>
  <td><img width=129 src="serial/TreeK-D=3-R=2-L=5-x-0.png"></img></td>
  <td><img width=129 src="serial/TreeK-D=3-R=2-L=6-x-0.png"></img></td>
  <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=7-x-0.png"></img></td> -->
  <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=8-x-0.png"></img></td> -->
  </tr>
  <tr>
  <th>project = Y</th>
  <td><img width=129 src="serial/TreeK-D=3-R=2-L=4-y-0.png"></img></td>
  <td><img width=129 src="serial/TreeK-D=3-R=2-L=5-y-0.png"></img></td>
  <td><img width=129 src="serial/TreeK-D=3-R=2-L=6-y-0.png"></img></td>
  <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=7-y-0.png"></img></td> -->
  <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=8-y-0.png"></img></td> -->
  </tr>
  <tr>
  <th>project = Z</th>
  <td><img width=129 src="serial/TreeK-D=3-R=2-L=4-z-0.png"></img></td>
  <td><img width=129 src="serial/TreeK-D=3-R=2-L=5-z-0.png"></img></td>
  <td><img width=129 src="serial/TreeK-D=3-R=2-L=6-z-0.png"></img></td>
  <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=7-z-0.png"></img></td> -->
  <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=8-z-0.png"></img></td> -->
  </tr>
  </table></br>


  <table>
  <tr>
  <th>coalesce = true</th>
  <th>levels = 4</th>
  <th>levels = 5</th>
  <th>levels = 6</th>
  <!-- <th>levels = 7</th> -->
  <!-- <th>levels = 8</th> -->
  </tr>
  <tr>
  <th>project = X</th>
  <td><img width=129 src="serial/TreeK-D=3-R=2-L=4-x-1.png"></img></td>
  <td><img width=129 src="serial/TreeK-D=3-R=2-L=5-x-1.png"></img></td>
  <td><img width=129 src="serial/TreeK-D=3-R=2-L=6-x-1.png"></img></td>
  <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=7-x-1.png"></img></td> -->
  <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=8-x-1.png"></img></td> -->
  </tr>
  <tr>
  <th>project = Y</th>
  <td><img width=129 src="serial/TreeK-D=3-R=2-L=4-y-1.png"></img></td>
  <td><img width=129 src="serial/TreeK-D=3-R=2-L=5-y-1.png"></img></td>
  <td><img width=129 src="serial/TreeK-D=3-R=2-L=6-y-1.png"></img></td>
  <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=7-y-1.png"></img></td> -->
  <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=8-y-1.png"></img></td> -->
  </tr>
  <tr>
  <th>project = Z</th>
  <td><img width=129 src="serial/TreeK-D=3-R=2-L=4-z-1.png"></img></td>
  <td><img width=129 src="serial/TreeK-D=3-R=2-L=5-z-1.png"></img></td>
  <td><img width=129 src="serial/TreeK-D=3-R=2-L=6-z-1.png"></img></td>
  <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=7-z-1.png"></img></td> -->
  <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=8-z-1.png"></img></td> -->
  </tr>
  </table></br>


  <h3>TreeK-D3-R4-L?</h3>

  <?php
  tests("Mesh","test_TreeK","test_TreeK-D3-R4-L4", "3D L=4 r=4");
  tests("Mesh","test_TreeK","test_TreeK-D3-R4-L6", "3D L=6 r=4");
  ?>

<table>
<tr>
<th></th>
<th colspan=2>coalesce = false</th>
  <th colspan=2>coalesce = true</th>
  </tr>
  <tr>
  <th></th>
  <th>levels = 4</th>
  <th>levels = 6</th>
  <!-- <th>levels = 8</th> -->
  <th>levels = 4</th>
  <th>levels = 6</th>
  <!-- <th>levels = 8</th> -->
  </tr>
  <tr>
  <th>project = X</th>
  <td><img width=129 src="serial/TreeK-D=3-R=4-L=4-x-0.png"></img></td>
  <td><img width=129 src="serial/TreeK-D=3-R=4-L=6-x-0.png"></img></td>
  <!-- <td><img width=129 src="serial/TreeK-D=3-R=4-L=8-x-0.png"></img></td> -->
  <td><img width=129 src="serial/TreeK-D=3-R=4-L=4-x-1.png"></img></td>
  <td><img width=129 src="serial/TreeK-D=3-R=4-L=6-x-1.png"></img></td>
  <!-- <td><img width=129 src="serial/TreeK-D=3-R=4-L=8-x-1.png"></img></td> -->
  </tr>
  <tr>
  <th>project =  Y</th>
  <td><img width=129 src="serial/TreeK-D=3-R=4-L=4-y-0.png"></img></td>
  <td><img width=129 src="serial/TreeK-D=3-R=4-L=6-y-0.png"></img></td>
  <!-- <td><img width=129 src="serial/TreeK-D=3-R=4-L=8-y-0.png"></img></td> -->
  <td><img width=129 src="serial/TreeK-D=3-R=4-L=4-y-1.png"></img></td>
  <td><img width=129 src="serial/TreeK-D=3-R=4-L=6-y-1.png"></img></td>
  <!-- <td><img width=129 src="serial/TreeK-D=3-R=4-L=8-y-1.png"></img></td> -->
  </tr>
  <tr>
  <th>project = Z</th>
  <td><img width=129 src="serial/TreeK-D=3-R=4-L=4-z-0.png"></img></td>
  <td><img width=129 src="serial/TreeK-D=3-R=4-L=6-z-0.png"></img></td>
  <!-- <td><img width=129 src="serial/TreeK-D=3-R=4-L=8-z-0.png"></img></td> -->
  <td><img width=129 src="serial/TreeK-D=3-R=4-L=4-z-1.png"></img></td>
  <td><img width=129 src="serial/TreeK-D=3-R=4-L=6-z-1.png"></img></td>
  <!-- <td><img width=129 src="serial/TreeK-D=3-R=4-L=8-z-1.png"></img></td> -->
  </br/>
  </body>
  </html>
