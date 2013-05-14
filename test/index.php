<html>
<head>
<title>Enzo-P / Cello Test Results</title>
<link href="cello.css" rel="stylesheet" type="text/css">
<?php
if (file_exists("COMPILING")) {
  echo "<meta http-equiv=\"refresh\" content=20>";
}
?>
   </head>
   <body>
   <h1>Enzo-P / Cello Test Results</h1>

   <h2><center><?php system ("hg branch") ?> branch </center> </h2>

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

function test_duration ($output_file) {
  if (file_exists($output_file)) {
    system("cat $output-file | awk '/UNIT TEST END/ {print $4};");
  } else {
    echo "<td></td>";
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

$types = array("charm","mpi");
  $num_types = sizeof($types);

  function tests($component,$testrun,$output,$test_name) {

  global $types;
  global $num_types;

  $input_file  = "../input/$output.in";
  $input_link  = "../input/$output.in";
  $input_html = "<a href=\"$input_link\">$output.in</a>";

  $source_file = "../src/$component/$testrun.cpp";
  $source_link = "../src/$component/$testrun.cpp";
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

function summary_incomplete_output ( $test_output, $executables)
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

    system("cat $output_files | awk 'BEGIN{e=0;}; /UNIT TEST END/ {e=e+1};END{if ($num_output_files==e) {print \"<td></td>\"} else {print \"<td class=fail>\"$num_output_files - e\"</td>\";}}'");

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
  summary_missing_output     ($test_output, $executables);
  summary_incomplete_output  ($test_output, $executables);
  summary_failed_tests       ($test_output, $executables);
  summary_unfinished_tests   ($test_output, $executables);
  summary_passed_tests       ($test_output, $executables);

  printf ("</tr>\n");
}

function test_table ($file_root,$size_array, $types)
{
  $show_flash = 0;
  echo "<table>";
  echo "<tr>";
  echo "<th>$file_root</th>";
  if ($show_flash) echo "<th>animation</th>";
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
	if ($show_flash) {
	  echo "<td>";
	  $swf_file = "$type/$file_root.swf"; 
	  $size_last = $size_array[sizeof($size_array)-1]; 
	  $png_file_last = "$type/$file_root-$size_last.png"; 
	  swf_movie($swf_file, 
		    $png_file_last, 
		    160); 
	  echo "</td>";
	}
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

function test_table_blocks ($file_root,$cycle_array, $types)
{
  echo "<table>";
  $rows = 2;
  $cols = 4;
  for ($i = 0; $i < sizeof($types); ++$i) {
     $type = $types[$i];
     echo "<tr><th>$file_root</th>";
     for ($index_cycle = 0; $index_cycle < sizeof($cycle_array); $index_cycle++) {
       $cycle = $cycle_array[$index_cycle];
       echo "<th colspan = $cols class=block>$cycle</th> ";
     }
     echo " </tr>";
     for ($row = 0; $row < $rows; $row++) {
       echo "<tr>";
       echo "<td>$type</td>";
       
       for ($index_cycle = 0; $index_cycle < sizeof($cycle_array); $index_cycle++) {
	 $cycle = $cycle_array[$index_cycle];
	 for ($col = 0; $col < $cols; $col++) {
	   $block = sprintf ("%08d-%08d-%08d",$rows - $row - 1,$col,0);
	   echo "<td class=block> <img src=$type/${file_root}-$cycle-block_$block.png width=80></img> </td>";
	 }
       }
       echo "</tr>";
     }
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
  
printf ( "<th colspan=$num_types class=fail>Missing</br>Executable</th>");
printf ("<th></th>");
printf ( "<th colspan=$num_types class=fail>Missing</br>Output</th>");
printf ("<th></th>");
printf ( "<th colspan=$num_types class=fail>Incomplete</br>Output</th>");
printf ("<th></th>");
printf ( "<th colspan=$num_types class=fail>Failed</br>Tests</th>");
printf ("<th></th>");
printf ( "<th colspan=$num_types class=unfinished>Unfinished</br>Tests</th>");
printf ("<th></th>");
printf ( "<th colspan=$num_types class=pass>Passed</br>Tests</th>");
printf ("<th></th>");
printf ( "</tr><tr>\n");

for ($k = 0; $k < 6; $k ++) {
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
	     array("method_ppm-1",
		   "method_ppm-8",
                   "checkpoint_ppm-1",
		   "checkpoint_ppm-8"),
	     array("enzo-p",  "enzo-p", "enzo-p", "enzo-p"));

test_summary("Enzo-PPML",
	     array("method_ppml-1",
		   "method_ppml-8"),
	     array("enzo-p",  "enzo-p"));

test_summary("Enzo-AMR", 
	     array("adapt-L1-P1", "adapt-L3-P1",),
	     array("enzo-p", "enzo-p"));

test_summary("Enzo-BC-2D", 
	     array("boundary_reflecting-2d",
		   "boundary_periodic-2d",
		   "boundary_outflow-2d"),
	     array("enzo-p", "enzo-p", "enzo-p", "enzo-p"));

test_summary("Enzo-BC-3D",
	     array("boundary_reflecting-3d",
		   "boundary_periodic-3d",
		   "boundary_outflow-3d"),
	     array("enzo-p", "enzo-p", "enzo-p", "enzo-p"));

test_summary("Enzo-IC", 
	     array("initial_png"),
	     array("enzo-p"));

test_summary("Enzo-output", 
	     array("output-stride-1","output-stride-2","output-stride-3"),
	     array("enzo-p","enzo-p","enzo-p"));

// Print row divider

printf ("<tr><th></th>");
for ($k = 0; $k < 6; $k ++) {
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
test_summary("Field",array("FieldBlock","FieldDescr","FieldFace","ItField"),
	     array("test_FieldBlock","test_FieldDescr","test_FieldFace","test_ItField")); 
test_summary("Io",array("ItReduce"),
	     array("test_ItReduce")); 
test_summary("Memory",array("Memory"),
	     array("test_Memory")); 
test_summary("Mesh",
	     array("Hierarchy",
		   "Block",
		   "Tree",
		   "TreeDensity",
		   "Node",
		   "NodeTrace",
		   "ItNode"),
	     array("test_Hierarchy",
		   "test_Block",
		   "test_Index",
		   "test_Tree",
		   "test_TreeDensity",
		   "test_Node",
		   "test_NodeTrace",
		   "test_ItNode")); 
test_summary("Monitor",array("Monitor"),
	     array("test_Monitor")); 
test_summary("Parallel",array("GroupProcess","Layout"),
	     array("test_GroupProcess","test_Layout")); 
test_summary("Parameters",array("Parameters"),
	     array("test_Parameters")); 
test_summary("Performance",array("Papi", "Performance","Timer"),
	     array("test_Papi","test_Performance","test_Timer")); 


printf ("</tr></table></br>\n");

//======================================================================

test_group("Enzo-PPM");

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


echo "<h3>PPM (serial) </h3>";

tests("Enzo","enzo-p","test_method_ppm-1","PPM 1 block");

test_table ("method_ppm-1",
	    array("000000","000200","000400"), $types);

//----------------------------------------------------------------------

echo "<h3>PPM (parallel) </h3>";

tests("Enzo","enzo-p","test_method_ppm-8","PPM 8 blocks");

test_table ("method_ppm-8",
	    array("000000","000200","000400"), $types);

//----------------------------------------------------------------------

echo "<h3>PPM (serial checkpoint / restart) </h3>";


tests("Enzo","enzo-p","test_checkpoint_ppm-1","test_restart_ppm-1","");

test_table ("restart_ppm-1",
	    array("000000","000200","000400"), $types);

//----------------------------------------------------------------------

echo "<h3>PPM (parallel checkpoint / restart) </h3>";


tests("Enzo","enzo-p","test_checkpoint_ppm-8","test_restart_ppm-8","");

test_table ("restart_ppm-8",
	    array("000000","000200","000400"), $types);

   //======================================================================


test_group("Enzo-PPML");

?>

Enzo-PPML tests serve to test basic PPML functionality in Enzo-P.  A
small high-density sphere is run for 50 cycles, first with
  one block (1,1,1) then eight blocks (2,2,2).

<?php

echo "<h3>PPML (serial) </h3>";

tests("Enzo","enzo-p","test_method_ppml-1","PPML 1 block");

test_table ("method_ppml-1-x",
	    array("0000","0010","0020","0030","0040"), $types);
test_table ("method_ppml-1-y",
	    array("0000","0010","0020","0030","0040"), $types);
test_table ("method_ppml-1-z",
	    array("0000","0010","0020","0030","0040"), $types);

echo "<h3>PPML (parallel) </h3>";

tests("Enzo","enzo-p","test_method_ppml-8","PPML 8 blocks");

test_table ("method_ppml-8-x",
	    array("0000","0010","0020","0030","0040"), $types);
test_table ("method_ppml-8-y",
	    array("0000","0010","0020","0030","0040"), $types);
test_table ("method_ppml-8-z",
	    array("0000","0010","0020","0030","0040"), $types);


//======================================================================

test_group("Enzo-AMR");

echo "<h3>2D Serial</h3>";

tests("Enzo","enzo-p","test_adapt-L1-P1","Level 1");

test_table ("adapt-L1-P1-mesh",
	    array("0000","0020","0040","0060","0080"), $types);
test_table ("adapt-L1-P1-density",
	    array("0000","0020","0040","0060","0080"), $types);

test_table ("adapt-L3-P1-mesh",
	    array("0000","0020","0040","0060","0080"), $types);
test_table ("adapt-L3-P1-density",
	    array("0000","0020","0040","0060","0080"), $types);


//======================================================================

test_group("Enzo-BC-2D");

echo "<h3>2D Reflecting</h3>";

tests("Enzo","enzo-p","test_boundary_reflecting-2d","Reflecting 2D");

test_table ("boundary_reflecting-2d",
	    array("0000","0100","0200","0300","0400"), $types);


//----------------------------------------------------------------------

echo "<h3>2D Periodic</h3>";

tests("Enzo","enzo-p","test_boundary_periodic-2d","Periodic 2D");

test_table ("boundary_periodic-2d",
	    array("0000","0100","0200","0300","0400"), $types);

//----------------------------------------------------------------------

echo "<h3>2D Outflow</h3>";
  
tests("Enzo","enzo-p","test_boundary_outflow-2d","Outflow 2D");

test_table ("boundary_outflow-2d",
	    array("0000","0100","0200","0300","0400"), $types);

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

echo "<h3>3D Outflow</h3>";

tests("Enzo","enzo-p","test_boundary_outflow-3d","Outflow 3D");

test_table ("boundary_outflow-3d",
	    array("0000","0020","0040","0060","0080"), $types);

//----------------------------------------------------------------------

test_group("Enzo-IC");


echo "<h3>png mask initial conditions</h3>";

tests("Enzo","enzo-p","test_initial_png","");

test_table ("initial_png",
	    array("00","10","20","30","40", "50"), $types);

//----------------------------------------------------------------------

test_group("Enzo-output");


echo "<h3>Process stride</h3>";

echo "<h3>Stride 1</h3>";

tests("Enzo","enzo-p","test_output-stride-1","");

test_table_blocks ("output-stride-1", array("00","10","20"),$types);

echo "<h3>Stride 2</h3>";

tests("Enzo","enzo-p","test_output-stride-2","");

test_table_blocks ("output-stride-2",  array("00","10","20"), $types);

echo "<h3>Stride 3</h3>";

tests("Enzo","enzo-p","test_output-stride-3","");

test_table_blocks ("output-stride-3",  array("00","10","20"), $types);

   //======================================================================

test_group("Disk");

tests("Cello","test_FileHdf5", "test_FileHdf5","");
tests("Cello","test_FileIfrit","test_FileIfrit","");

   //----------------------------------------------------------------------

test_group("Error");

tests("Cello","test_Error","test_Error","");


   //----------------------------------------------------------------------

test_group("Field");

tests("Cello","test_FieldDescr","test_FieldDescr","");
tests("Cello","test_FieldBlock","test_FieldBlock","");
tests("Cello","test_FieldFace","test_FieldFace","");
tests("Cello","test_ItField","test_ItField","");

   //----------------------------------------------------------------------

test_group("Io");

tests("Cello","test_ItReduce", "test_ItReduce","");

   //----------------------------------------------------------------------

test_group("Memory");

tests("Cello","test_Memory","test_Memory","");


   //----------------------------------------------------------------------

test_group("Mesh");

tests("Cello","test_Hierarchy","test_Hierarchy",""); 
tests("Cello","test_Block","test_Block",""); 
tests("Cello","test_Index","test_Index",""); 
tests("Cello","test_Tree","test_Tree",""); 

printf ("<img width=257 src=\"mpi/test_tree_1-initial.png\"></img>\n");
printf ("<img width=257 src=\"mpi/test_tree_2-balanced.png\"></img>\n");
printf ("<img width=257 src=\"mpi/test_tree_3-merged.png\"></img></br>\n");

tests("Cello","test_TreeDensity","test_TreeDensity",""); 

printf ("<img width=257 src=\"mpi/density_xy_1-initial.png\"></img>\n");
printf ("<img width=257 src=\"mpi/density_xy_2-balanced.png\"></img>\n");
printf ("<img width=257 src=\"mpi/density_xy_3-coalesced.png\"></img></br>\n");

printf ("<img width=257 src=\"mpi/density_3d_1-initial.png\"></img>\n");
printf ("<img width=257 src=\"mpi/density_3d_2-balanced.png\"></img>\n");
printf ("<img width=257 src=\"mpi/density_3d_3-coalesced.png\"></img></br>\n");

tests("Cello","test_Node","test_Node",""); 
tests("Cello","test_NodeTrace","test_NodeTrace",""); 
tests("Cello","test_ItNode","test_ItNode",""); 

   //----------------------------------------------------------------------

test_group("Monitor");

tests("Cello","test_Monitor","test_Monitor","");

// printf ("<img src=\"monitor_image_1.png\"></img>\n");
// printf ("<img src=\"monitor_image_2.png\"></img>\n");
// printf ("<img src=\"monitor_image_3.png\"></img>\n");
// printf ("<img src=\"monitor_image_4.png\"></img>\n");

   //----------------------------------------------------------------------

test_group("Parallel");

tests("Cello","test_GroupProcess","test_GroupProcess","");
tests("Cello","test_Layout","test_Layout","");

   //----------------------------------------------------------------------

test_group("Parameters");

tests("Cello","test_Parameters","test_Parameters","");

   //----------------------------------------------------------------------

test_group("Performance");

tests("Cello","test_Performance","test_Performance","");
tests("Cello","test_Papi",       "test_Papi","");
tests("Cello","test_Timer",       "test_Timer","");

/* <hr> */
/* <h2>Mesh Tests (Prototype Code)</h2> */

/*   <h3>TreeK-D2-R2-L?</h3> */

/*   <?php */
/*   tests("Cello","test_TreeK","test_TreeK-D2-R2-L6", "2D L=6 r=2"); */
/*   tests("Cello","test_TreeK","test_TreeK-D2-R2-L7", "2D L=7 r=2"); */
/*   tests("Cello","test_TreeK","test_TreeK-D2-R2-L8", "2D L=8 r=2"); */
/*   ?> */

/* <table> */
/* <tr> */
/* <th>coalesce</th> */
/* <th>levels = 6</th> */
/*   <th>levels = 7</th> */
/*   <th>levels = 8</th> */
/*   </tr> */
/*   <tr> */
/*   <th>false</th> */
/*   <td><img width=257 src="serial/TreeK-D=2-R=2-L=6-0.png"></img></td> */
/*   <td><img width=257 src="serial/TreeK-D=2-R=2-L=7-0.png"></img></td> */
/*   <td><img width=257 src="serial/TreeK-D=2-R=2-L=8-0.png"></img></td> */
/*   </tr> */
/*   <tr> */
/*   <th>true</th> */
/*   <td><img width=257 src="serial/TreeK-D=2-R=2-L=6-1.png"></img></td> */
/*   <td><img width=257 src="serial/TreeK-D=2-R=2-L=7-1.png"></img></td> */
/*   <td><img width=257 src="serial/TreeK-D=2-R=2-L=8-1.png"></img></td> */
/*   </tr> */
/*   </table></br> */

/*   <h3>TreeK-D2-R4-L?</h3> */

/*   <?php */
/*   tests("Cello","test_TreeK","test_TreeK-D2-R4-L6", "2D L=6 r=4"); */
/*   tests("Cello","test_TreeK","test_TreeK-D2-R4-L8", "2D L=8 r=4"); */
/*   ?> */

/* <table> */
/* <tr> */
/* <th>coalesce</th> */
/* <th>levels = 6</th> */
/*   <th>levels = 8</th> */
/*   <!-- <th>levels = 10</th> --> */
/*   </tr> */
/*   <tr> */
/*   <th>false</th> */
/*   <td><img width=257 src="serial/TreeK-D=2-R=4-L=6-0.png"></img></td> */
/*   <td><img width=257 src="serial/TreeK-D=2-R=4-L=8-0.png"></img></td> */
/*   <!-- <td><img width=257 src="serial/TreeK-D=2-R=4-L=10-0.png"></img></td> --> */
/*   </tr> */
/*   <tr> */
/*   <th>true</th> */
/*   <td><img width=257 src="serial/TreeK-D=2-R=4-L=6-1.png"></img></td> */
/*   <td><img width=257 src="serial/TreeK-D=2-R=4-L=8-1.png"></img></td> */
/*   <!-- <td><img width=257 src="serial/TreeK-D=2-R=4-L=10-1.png"></img></td> --> */
/*   </tr> */
/*   </table></br> */

/*   <h3>TreeK-D3-R2-L?</h3> */

/*   <?php */
/*   tests("Cello","test_TreeK","test_TreeK-D3-R2-L4", "3D L=4 r=2"); */
/*   tests("Cello","test_TreeK","test_TreeK-D3-R2-L5", "3D L=5 r=2"); */
/*   tests("Cello","test_TreeK","test_TreeK-D3-R2-L6", "3D L=6 r=2"); */
/*    ?> */

/* <table> */
/* <tr> */
/* <th>coalesce = false</th> */
/*   <th>levels = 4</th> */
/*   <th>levels = 5</th> */
/*   <th>levels = 6</th> */
/*   <!-- <th>levels = 7</th> --> */
/*   <!-- <th>levels = 8</th> --> */
/*   </tr> */
/*   <tr> */
/*   <th>project = X</th> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=2-L=4-x-0.png"></img></td> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=2-L=5-x-0.png"></img></td> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=2-L=6-x-0.png"></img></td> */
/*   <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=7-x-0.png"></img></td> --> */
/*   <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=8-x-0.png"></img></td> --> */
/*   </tr> */
/*   <tr> */
/*   <th>project = Y</th> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=2-L=4-y-0.png"></img></td> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=2-L=5-y-0.png"></img></td> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=2-L=6-y-0.png"></img></td> */
/*   <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=7-y-0.png"></img></td> --> */
/*   <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=8-y-0.png"></img></td> --> */
/*   </tr> */
/*   <tr> */
/*   <th>project = Z</th> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=2-L=4-z-0.png"></img></td> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=2-L=5-z-0.png"></img></td> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=2-L=6-z-0.png"></img></td> */
/*   <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=7-z-0.png"></img></td> --> */
/*   <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=8-z-0.png"></img></td> --> */
/*   </tr> */
/*   </table></br> */


/*   <table> */
/*   <tr> */
/*   <th>coalesce = true</th> */
/*   <th>levels = 4</th> */
/*   <th>levels = 5</th> */
/*   <th>levels = 6</th> */
/*   <!-- <th>levels = 7</th> --> */
/*   <!-- <th>levels = 8</th> --> */
/*   </tr> */
/*   <tr> */
/*   <th>project = X</th> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=2-L=4-x-1.png"></img></td> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=2-L=5-x-1.png"></img></td> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=2-L=6-x-1.png"></img></td> */
/*   <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=7-x-1.png"></img></td> --> */
/*   <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=8-x-1.png"></img></td> --> */
/*   </tr> */
/*   <tr> */
/*   <th>project = Y</th> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=2-L=4-y-1.png"></img></td> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=2-L=5-y-1.png"></img></td> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=2-L=6-y-1.png"></img></td> */
/*   <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=7-y-1.png"></img></td> --> */
/*   <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=8-y-1.png"></img></td> --> */
/*   </tr> */
/*   <tr> */
/*   <th>project = Z</th> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=2-L=4-z-1.png"></img></td> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=2-L=5-z-1.png"></img></td> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=2-L=6-z-1.png"></img></td> */
/*   <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=7-z-1.png"></img></td> --> */
/*   <!-- <td><img width=129 src="serial/TreeK-D=3-R=2-L=8-z-1.png"></img></td> --> */
/*   </tr> */
/*   </table></br> */


/*   <h3>TreeK-D3-R4-L?</h3> */

/*   <?php */
/*   tests("Cello","test_TreeK","test_TreeK-D3-R4-L4", "3D L=4 r=4"); */
/*   tests("Cello","test_TreeK","test_TreeK-D3-R4-L6", "3D L=6 r=4"); */
/*   ?> */

/* <table> */
/* <tr> */
/* <th></th> */
/* <th colspan=2>coalesce = false</th> */
/*   <th colspan=2>coalesce = true</th> */
/*   </tr> */
/*   <tr> */
/*   <th></th> */
/*   <th>levels = 4</th> */
/*   <th>levels = 6</th> */
/*   <!-- <th>levels = 8</th> --> */
/*   <th>levels = 4</th> */
/*   <th>levels = 6</th> */
/*   <!-- <th>levels = 8</th> --> */
/*   </tr> */
/*   <tr> */
/*   <th>project = X</th> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=4-L=4-x-0.png"></img></td> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=4-L=6-x-0.png"></img></td> */
/*   <!-- <td><img width=129 src="serial/TreeK-D=3-R=4-L=8-x-0.png"></img></td> --> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=4-L=4-x-1.png"></img></td> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=4-L=6-x-1.png"></img></td> */
/*   <!-- <td><img width=129 src="serial/TreeK-D=3-R=4-L=8-x-1.png"></img></td> --> */
/*   </tr> */
/*   <tr> */
/*   <th>project =  Y</th> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=4-L=4-y-0.png"></img></td> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=4-L=6-y-0.png"></img></td> */
/*   <!-- <td><img width=129 src="serial/TreeK-D=3-R=4-L=8-y-0.png"></img></td> --> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=4-L=4-y-1.png"></img></td> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=4-L=6-y-1.png"></img></td> */
/*   <!-- <td><img width=129 src="serial/TreeK-D=3-R=4-L=8-y-1.png"></img></td> --> */
/*   </tr> */
/*   <tr> */
/*   <th>project = Z</th> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=4-L=4-z-0.png"></img></td> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=4-L=6-z-0.png"></img></td> */
/*   <!-- <td><img width=129 src="serial/TreeK-D=3-R=4-L=8-z-0.png"></img></td> --> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=4-L=4-z-1.png"></img></td> */
/*   <td><img width=129 src="serial/TreeK-D=3-R=4-L=6-z-1.png"></img></td> */
/*   <!-- <td><img width=129 src="serial/TreeK-D=3-R=4-L=8-z-1.png"></img></td> --> */
?>
  </br/>
  </body>
  </html>
