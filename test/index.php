<html>
<head>
<title>Enzo-P / Cello Test Results</title>
<link href="cello.css" rel="stylesheet" type="text/css">

   <?php
   if (file_exists("STATUS")) {
     echo "<meta http-equiv=\"refresh\" content=20>";
   }
   ?>
   </head>
   <body>
   <h1>Enzo-P / Cello Test Results</h1>

   <h3><?php system ("hg branch") ?> branch </br> revision <?php system ("hg id -n") ?> </h3>
 <?php

     //----------------------------------------------------------------------

   function test_group($testgroup) {
     printf ("<hr><a name=\"$testgroup\"><h2>$testgroup</h2>");
     }

     //----------------------------------------------------------------------

   function test_input ($input_file) {
     $file="$input_file.in";
     if (file_exists($file)) {
       $input_html = "<a href=\"$file\">$input.unit.in</a>";
       printf ("<td class=pass><a href='$file'>input</a></td>");
     } else {
       echo "<td></td>";
     }
   }

//----------------------------------------------------------------------

   function test_output ($output_file) {
     $file="$output_file";
     if (file_exists($file)) {
       $output_html = "<a href=\"$file\">$output.unit</a>";
       system("cat $file | awk 'BEGIN{c=0}; /END CELLO/ {c=1}; END{ if (c==0) print \"<td class=fail><a href='$file'>output</a></td>\"; if (c!=0) print \"<td class=pass><a href='$file'>output</a></td>\"}'");
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
    echo "<td class=pass>";
    system("awk '/UNIT TEST END/ {print \"scale=0; \" $4 \"/1\"}' < $output_file | bc"); echo "s";
    echo "</td>";
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

$types = array("charm");
$num_types = sizeof($types);

function tests($component,$testrun,$output,$test_name,$dir) {

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
    echo "   <th>Input</th>";
    echo "   <th>Output</th>";
    echo "   <th>Date</th>";
    echo "   <th>Time</th>";
    echo "   <th>Duration</th>";
    echo "   <th>Failed</th>";
    echo "   <th>Passed</th>";
    //    echo "   <th>Unfinished</th>";
    echo "</tr>\n";

    //--------------------------------------------------

    for ($i = 0; $i<$num_types; ++$i) {

      echo "<tr>\n";

      //      echo "<th> $types[$i] </th>";

      $output_file = "$output.unit";

      if ($dir != "") {
	$file = "$dir/$output_file";
      } else {
	$file = "$output_file";
      }
      test_input      ($file);
      test_output     ($file);
      test_date       ($file);
      test_time       ($file);
      test_duration   ($file);
      test_failed     ($file);
      test_passed     ($file);
      //      test_unfinished ($file);

      echo "</tr>\n";
    }
    echo "</tr></table></br/>\n";

    echo "<table><tr>";
    for ($i = 0; $i<$num_types; ++$i) {
      $type = $types[$i];
      $output_file = "../test/$dir/$output.unit";
      if (file_exists($output_file)) {
	test($type,$output_file,"FAIL");
      }
    }
    echo "</tr></table>";
  }
};

//----------------------------------------------------------------------

function test($type,$output,$type) {
  $ltype = strtolower($type);

  $cols = "$4,$6,$7,$8,$9,$10";
  $rowtext = "</tr><tr>";

  $count = exec("cat $output | grep $type | grep '0/' | wc -l");
  if ($count != 0) {
    echo "<th class=$type><strong>${ltype}ed</strong></th> ";
    system ("grep '0/' $output | sort | uniq | awk 'BEGIN {c=1}; / $type /{split($3,a,\"\/\"); print \"<td class=$type> \",$cols , \" </td>\"; c=c+1}; {if (c==5) {c=0; print \"$rowtext\"}}'");
    echo "</tr><tr></tr>";
  }
};

?>

<hr>

<?php

 //----------------------------------------------------------------------

function summary_missing_executable ($test_output, $executables, $state, $dir, $valid)
{
  global $types;
  global $num_types;

  for ($i = 0; $i<$num_types; ++$i) {

    $type = $types[$i];

    $count_missing = 0;

    for ($test = 0; $test<sizeof($executables); ++$test) {
      $exe = $executables[$test];
      $bin = "../bin/$exe";
      if (! file_exists($bin)) {
	++ $count_missing ;
      }
    }

    if ($count_missing == 0 || ! $valid) { 
      printf ("<td></td>");
    } else {
      printf ("<td class=noexec>$count_missing</td>");
      $valid = false;
    }
  }
  return $valid;
}

//----------------------------------------------------------------------

function summary_missing_output ($test_output, $executables, $state, $dir,$valid)
{
  global $types;
  global $num_types;

  for ($i = 0; $i<$num_types; ++$i) {

    $count_missing = 0;
    for ($test = 0; $test<sizeof($test_output); ++$test) {
      $output = "../$dir/test_$test_output[$test].unit";
      if (! file_exists($output)) {
	++ $count_missing;
      }
    }
    if ($count_missing == 0 || ! $valid) {
      printf ("<td></td>");
    } else {
      printf ("<td class=noout>$count_missing</td>");
      $valid = false;
    }
  }
  return $valid;
}

//----------------------------------------------------------------------

function summary_incomplete_output ( $test_output, $executables, $state, $dir,$valid)
{
  global $types;
  global $num_types;

  for ($i = 0; $i<$num_types; ++$i) {

    $output_files = "";
    $num_output_files = 0;
    for ($test = 0; $test<sizeof($test_output); ++$test) {
      $output = "../$dir/test_$test_output[$test].unit";
      $output_files = "$output_files $output";
      ++$num_output_files;
    }

    /* if (! $valid) { */
    /*   printf ("<td></td>"); */
    /* } else { */
    passthru("cat $output_files | awk 'BEGIN{e=0;}; /END CELLO/ {e=e+1};END{if ($num_output_files==e) {print \"<td></td>\"} else {print \"<td class=incomplete>\"$num_output_files - e\"</td>\";}}'");
    /* } */

  }
  return $valid;
}

//----------------------------------------------------------------------

function summary_failed_tests ($test_output, $executables, $state, $dir,$valid)
{
  global $types;
  global $num_types;

  for ($i = 0; $i<$num_types; ++$i) {

    $output_files = "";
    for ($test = 0; $test<sizeof($test_output); ++$test) {
      $output = "../$dir/test_$test_output[$test].unit";
      $output_files = "$output_files $output";
    }
    if (! $valid) {
      printf ("<td></td>");
    } else {
      system("grep '0/' $output_files | awk 'BEGIN {c=0}; /FAIL/{c=c+1}; END{if (c==0) {print \"<td></td>\"} else {print \"<td class=fail>\",c,\"</td>\";}} '");
    }
  }
  return $valid;
}

//----------------------------------------------------------------------

function summary_unfinished_tests ($test_output, $executables, $state, $dir,$valid)
{
  global $types;
  global $num_types;

  for ($i = 0; $i<$num_types; ++$i) {

    $output_files = "";
    for ($test = 0; $test<sizeof($test_output); ++$test) {
      $output = "../$dir/test_$test_output[$test].unit";
      $output_files = "$output_files $output";
    }

    /* if (! $valid) { */
    /*   printf ("<td></td>"); */
    /* } else { */
      system("grep '0/' $output_files | awk 'BEGIN {c=0}; /incomplete/{c=c+1}; END{if (c==0) {print \"<td></td>\"} else {print \"<td class=unfinished>\",c,\"</td>\";}} '");
    /* } */
  }
  return $valid;
}

//----------------------------------------------------------------------

function summary_passed_tests ($test_output, $executables, $state, $dir,$valid)
{

  global $types;
  global $num_types;

  for ($i = 0; $i<$num_types; ++$i) {

    $output_files = "";
    for ($test = 0; $test<sizeof($test_output); ++$test) {
      $output = $test_output[$test];
      $output_files = "$output_files ../$dir/test_$output.unit";
    }
    system("grep '0/' $output_files | awk 'BEGIN {c=0}; /pass/{c=c+1}; END{if (c==0) {print \"<td></td>\"} else {print \"<td class=pass>\",c,\"</td>\";}} '");

  }
  return $valid;
}

//----------------------------------------------------------------------

function test_summary($component,$test_output,$executables, $dir)
{

  printf ("<tr><th><a href=\"#$component\">$component</a></th>\n");

  global $types;

  // Missing executable

  $state = "running";

  $valid = true;
  $valid = summary_missing_executable 
    ($test_output, $executables, $state, $dir,$valid);
  $valid = summary_missing_output
    ($test_output, $executables, $state, $dir,$valid);
  $valid = summary_incomplete_output
    ($test_output, $executables, $state, $dir,$valid);
  $valid = summary_failed_tests
    ($test_output, $executables, $state, $dir,$valid);
  $valid = summary_passed_tests
    ($test_output, $executables, $state, $dir,$valid);

  /* printf ("<th class=divider></th>"); */
  /* $valid = summary_unfinished_tests */
  /*   ($test_output, $executables, $state, $dir,$valid); */

  printf ("</tr>\n");
}


function begin_hidden ($id, $title)
{  
  echo "<a href=\"#\" onclick=\"document.getElementById('$id').style.display='block'; return false;\"><h3>$title</h3></a></p> ";
  echo "<div style=\"display: none;\" id=\"$id\">";
  echo "<a href=\"#\" onclick=\"document.getElementById('$id').style.display='none'; return false;\">[hide]</a></p>";
}

function end_hidden ($id)
{
  echo "<a href=\"#\" onclick=\"document.getElementById('$id').style.display='none'; return false;\">[hide]";
  echo "</div>";
}

function test_table ($file_root,$size_array, $types)
{
  $show_flash = 0;

  echo "<table>";
  echo "<tr>";
  //  echo "<th>$file_root</th>";
  if ($show_flash) echo "<th>animation</th>";
  for ($j = 0; $j < sizeof($size_array); ++$j) {
    $size = $size_array[$j];
    printf ("<th>$size</th>\n");
  }
  echo "</tr>";
  for ($i = 0; $i < sizeof($types); ++$i) {
    echo "<tr>";
    $type = $types[$i];
    //     	printf ("<th>$type</th>\n"); 
    // Show movie file if available
    if ($show_flash) {
      echo "<td>";
      $swf_file = "$type/$file_root.swf"; 
      $size_last = $size_array[sizeof($size_array)-1]; 
      $png_file_last = "$file_root-$size_last.png"; 
      swf_movie($swf_file, 
		$png_file_last, 
		160); 
      echo "</td>";
    }
    // Show available image frames
    for ($j = 0; $j < sizeof($size_array); ++$j) {
      $size = $size_array[$j];
      $png_file = "$file_root-$size.png"; 

      printf ("<td><img width=160 src=$png_file></img></td>\n");  
    }  
    echo "</tr>";  
  }
  echo "</table></br>";
}

function binary ($value,$count)
{
  $strval = "";
  while ($count > 0) {
    $bit = (int)$value % 2;
    $strval = "$bit".$strval;
    $value = $value >> 1;
    $count = $count - 1;
  }
  return $strval; 
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
	  $rbin = binary($rows - $row - 1,2);
	  $cbin = binary($col,2);
	  echo "<td class=block> <img src=${file_root}-$cycle-B${rbin}_${cbin}.png width=80></img> </td>";
	}
      }
      echo "</tr>";
    }
  }
  echo "</table></br>";
}

function row_divider($num_types)
{
  printf ("<tr>");
  printf ("<th class=divider></th>");
  for ($k = 0; $k < 5; $k ++) {
    printf ("<th class=divider></th>");
    for ($i = 0; $i < $num_types; ++$i) {
      printf ("<th class=divider> </th>");
    }
  }
  printf ("<th class=divider></th>");
  printf ("</tr>\n");
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

// Print the top row of the summary table

if (file_exists("STATUS"))  {
  printf ( "<th class=compiling ><strong>Running");
  printf ("</strong></th>");
} else { 
  printf ("<th></th>\n"); 
}
printf ( "<th colspan=$num_types class=fail>No program</th>");
printf ( "<th colspan=$num_types class=fail>No output</th>");
printf ( "<th colspan=$num_types class=fail>Some output</th>");
printf ( "<th colspan=$num_types class=fail>Failed </th>");
printf ( "<th colspan=$num_types class=pass>Passed </th>");
/* printf ("<th class=divider></th>"); */
/* printf ( "<th colspan=$num_types class=unfinished>Unfinished</br>Tests</th>"); */
printf ( "</tr>\n");

//----------------------------------------------------------------------
//----------------------------------------------------------------------

if (file_exists("STATUS"))  {
  printf ("<tr><th class=compiling>");
  system ("awk '{i=index($1,\"test_\"); print substr($1,i+5,length($1)-i-9)}' STATUS");
}  else {
  printf ("<tr><th>\n");
}

printf ("</th><td class=center colspan=5><em><a href=\"#enzop\">Enzo-P application tests</a></em></td></tr>\n");

test_summary("Method: ppm",
	     array("method_ppm-1","method_ppm-8"),
	     array("enzo-p",  "enzo-p"),'test');

test_summary("Method: ppml",
	     array("method_ppml-1","method_ppml-8"),
	     array("enzo-p",  "enzo-p"),'test');

test_summary("Method: heat",
	     array("method_heat-1","method_heat-8"),
	     array("enzo-p",  "enzo-p"),'test');

test_summary("Method: gravity",
	     array("method_gravity_cg-1","method_gravity_cg-8"),
	     array("enzo-p",  "enzo-p"),'test');

test_summary("Checkpoint",
	     array("checkpoint_ppm-1","checkpoint_ppm-8"),
	     array("enzo-p",  "enzo-p"),'test');

test_summary("Mesh", 
	     array("mesh-balanced"),
	     array("enzo-p"),'test');

test_summary("Balance", 
	     array("balance_none",
		   "balance_rand_cent",
		   "balance_greedy",
		   "balance_refine",
		   "balance_rotate",
		   "balance_hybrid"),
	     array("enzo-p", "enzo-p", "enzo-p", "enzo-p",
		   "enzo-p", "enzo-p"),'test/Balance');

test_summary("Boundary", 
	     array("boundary_reflecting-2d",
		   "boundary_periodic-2d",
		   "boundary_outflow-2d",
		   "boundary_reflecting-3d",
		   "boundary_periodic-3d",
		   "boundary_outflow-3d"),
	     array("enzo-p", "enzo-p", "enzo-p",
		   "enzo-p", "enzo-p", "enzo-p"),'test');

test_summary("Initial", 
	     array("initial_png"),
	     array("enzo-p"),'test');

test_summary("Output", 
	     array("output-stride-1","output-stride-2","output-stride-4"),
	     array("enzo-p","enzo-p","enzo-p"),'test');

//----------------------------------------------------------------------
// row_divider($num_types);
printf ("<tr><td class=center></td><td class=center colspan=5><em><a href=\"#cello\">Cello unit tests</a></em></td></tr>\n");
//----------------------------------------------------------------------


test_summary("Disk",array(     "FileHdf5",     "FileIfrit"),
	     array("test_FileHdf5","test_FileIfrit"),'test'); 
test_summary("Error",array(    "Error"),
	     array("test_Error"),'test'); 
test_summary("Field",
	     array(     "Field",      "FieldData",     "FieldDescr",     "FieldFace",     "ItField",      "Grouping"),
	     array("test_Field", "test_FieldData","test_FieldDescr","test_FieldFace","test_ItField", "test_Grouping"),
	     'test'); 
test_summary("Memory",array("Memory"),
	     array("test_Memory"),'test'); 
test_summary("Mesh",
	     array("Hierarchy",
		   "Data",
		   "Index",
		   "Tree",
		   "ItFace",
		   "TreeDensity",
		   "Node",
		   "NodeTrace",
		   "ItNode"),
	     array("test_Hierarchy",
		   "test_Data",
		   "test_Index",
		   "test_Tree",
		   "test_ItFace",
		   "test_TreeDensity",
		   "test_Node",
		   "test_NodeTrace",
		   "test_ItNode"),'test'); 
test_summary("Monitor",array("Monitor"),
	     array("test_Monitor"),'test'); 
test_summary("Parameters",
	     array("Parameters"),
	     array("test_Parameters"),'test'); 
test_summary("Particle",array("Particle"),
	     array("test_Particle"),'test'); 
test_summary("Performance",array("Papi", "Performance","Timer"),
	     array("test_Papi","test_Performance","test_Timer"),'test'); 
test_summary("Problem",array("Mask","Refresh","Value"),
	     array("test_Mask","test_Refresh","test_Value"),'test'); 
test_summary("Schedule",array("Schedule"),
	     array("test_Schedule"),'test'); 
test_summary("Colormap",array("Colormap"),
	     array("test_Colormap"),'test'); 

printf ("</tr></table></br>\n");

   /* <code>Start: <?php system ("cat START") ?> </code><br> */
   /* <code>Stop:&nbsp; <?php system ("cat STOP") ?></code> */

//======================================================================
     
     
     echo "<a name=\"enzop\"><h1>Enzo-P application tests</h1>";

test_group("Method: ppm");

?>

Method-PPM tests serve to test basic PPM functionality in Enzo-P.  A
small implosion problem is run for 400 cycles, first with
  one block (1,1) then eight blocks (2,4).
<p>

  Currently, "serial" results are incorrect for multiple blocks, which 
  is to be expected.  There are slight discrepencies in parallel runs
  eight blocks because the final time after 400 cycles does not exactly 
  match the time for the serial runs.  The results look qualitatively 
  correct however, even at time 2.5 for 400<sup>2</sup> (over 13000 cycles). 
  </p>

  <?php



//----------------------------------------------------------------------

begin_hidden("method_ppm-1", "PPM (serial)");

tests("Enzo","enzo-p","test_method_ppm-1","PPM 1 block","");

test_table ("method_ppm-1",
	    array("000000","000200","000400"), $types);

end_hidden("method_ppm-1");

//----------------------------------------------------------------------

begin_hidden("method_ppm-8", "PPM (parallel)");

tests("Enzo","enzo-p","test_method_ppm-8","PPM 8 blocks","");

?>
See <a href="http://client64-249.sdsc.edu/cello-bug/show_bug.cgi?id=19">Bug #19</a> for "final time" discrepency between serial and parallel PPM runs. </p>
<?php

test_table ("method_ppm-8",
	    array("000000","000200","000400"), $types);

end_hidden("method_ppm-8");

//======================================================================


test_group("Method: ppml");

?>

Method-PPML tests serve to test basic PPML functionality in Enzo-P.  A
small high-density sphere is run for 50 cycles, first with one 
block (1,1,1) then eight blocks (2,2,2).

  <?php

  begin_hidden("method_ppml-1", "PPML (serial)");

tests("Enzo","enzo-p","test_method_ppml-1","PPML 1 block","");

test_table ("method_ppml-1-x",
	    array("0000","0010","0020","0030","0040"), $types);
test_table ("method_ppml-1-y",
	    array("0000","0010","0020","0030","0040"), $types);
test_table ("method_ppml-1-z",
	    array("0000","0010","0020","0030","0040"), $types);

end_hidden ("method_ppml-1");

begin_hidden("method_ppml-8", "PPML (parallel)");

tests("Enzo","enzo-p","test_method_ppml-8","PPML 8 blocks","");

test_table ("method_ppml-8-x",
	    array("0000","0010","0020","0030","0040"), $types);
test_table ("method_ppml-8-y",
	    array("0000","0010","0020","0030","0040"), $types);
test_table ("method_ppml-8-z",
	    array("0000","0010","0020","0030","0040"), $types);

end_hidden ("method_ppml-8");

//======================================================================

test_group("Method: heat");

?>

Method-heat tests serve to test basic functionality of the "heat" method
in Enzo-P.

</p>

<?php


  begin_hidden("method_heat-1", "HEAT (serial)");

tests("Enzo","enzo-p","test_method_heat-1","HEAT 1 block","");

test_table ("method_heat-temp-1",
	    array("000000","000200","000400"), $types);
test_table ("method_heat-mesh-1",
	    array("000000","000200","000400"), $types);

end_hidden ("method_heat-1");

  begin_hidden("method_heat-8", "HEAT (parallel)");

tests("Enzo","enzo-p","test_method_heat-8","HEAT 8 block","");

test_table ("method_heat-temp-8",
	    array("000000","000200","000400"), $types);
test_table ("method_heat-mesh-8",
	    array("000000","000200","000400"), $types);

end_hidden ("method_heat-8");

//======================================================================

test_group("Method: gravity");

?>

Method-gravity tests serve to test basic functionality of the "gravity_cg" method
in Enzo-P.

</p>

<?php


  begin_hidden("method_gravity_cg-1", "GRAVITY (serial)");

tests("Enzo","enzo-p","test_method_gravity_cg-1","GRAVITY_CG 1 block","");

test_table ("method_gravity_cg-1",
	    array("mesh-000000","mesh-000010","mesh-000020","mesh-000030","mesh-000040","mesh-000050"), $types);
test_table ("method_gravity_cg-1",
	    array("rho-000000","rho-000010","rho-000020","rho-000030","rho-000040","rho-000050"), $types);
test_table ("method_gravity_cg-1",
	    array("phi-000000","phi-000010","phi-000020","phi-000030","phi-000040","phi-000050"), $types);
test_table ("method_gravity_cg-1",
	    array("ax-000000","ax-000010","ax-000020","ax-000030","ax-000040","ax-000050"), $types);
test_table ("method_gravity_cg-1",
	    array("ay-000000","ay-000010","ay-000020","ay-000030","ay-000040","ay-000050"), $types);

end_hidden("method_gravity_cg-1");

  begin_hidden("method_gravity_cg-8", "GRAVITY (parallel)");

tests("Enzo","enzo-p","test_method_gravity_cg-8","GRAVITY_CG 8 block","");

test_table ("method_gravity_cg-8",
	    array("mesh-000000","mesh-000010","mesh-000020","mesh-000030","mesh-000040","mesh-000050"), $types);
test_table ("method_gravity_cg-8",
	    array("rho-000000","rho-000010","rho-000020","rho-000030","rho-000040","rho-000050"), $types);
test_table ("method_gravity_cg-8",
	    array("phi-000000","phi-000010","phi-000020","phi-000030","phi-000040","phi-000050"), $types);
test_table ("method_gravity_cg-8",
	    array("ax-000000","ax-000010","ax-000020","ax-000030","ax-000040","ax-000050"), $types);
test_table ("method_gravity_cg-8",
	    array("ay-000000","ay-000010","ay-000020","ay-000030","ay-000040","ay-000050"), $types);

end_hidden("method_gravity_cg-8");

//======================================================================

test_group("Checkpoint");

begin_hidden("checkpoint_ppm-1","Checkpoint/Restart (serial)");

tests("Enzo","enzo-p","test_checkpoint_ppm-1","Checkpoint P=1","");
tests("Enzo","enzo-p","test_restart_ppm-1","Restart P=1","");
test_table ("checkpoint_ppm-1",  array("000010","000020"), $types);

end_hidden("checkpoint_ppm-1");

//----------------------------------------------------------------------

begin_hidden("checkpoint_ppm-8","Checkpoint/Restart (parallel)");


tests("Enzo","enzo-p","test_checkpoint_ppm-8","Checkpoint P=8","");
tests("Enzo","enzo-p","test_restart_ppm-8","Restart P=8","");
test_table ("checkpoint_ppm-1",  array("000010","000020"), $types);

end_hidden("checkpoint_ppm-8");

//======================================================================

test_group("Mesh");

begin_hidden ("mesh-balanced", "Adapt (serial)");

tests("Enzo","enzo-p","test_mesh-balanced","balanced","");

test_table ("mesh-balanced", array("mesh.000","de.000","te.000","vx.000","vy.000"), $types);
test_table ("mesh-balanced", array("mesh.100","de.100","te.100","vx.100","vy.100"), $types);

end_hidden ("mesh-balanced");

/* //====================================================================== */

/* test_group("Enzo-AMR"); */

begin_hidden("adapt_L5", "Adapt (parallel)");

tests("Enzo","enzo-p","test_adapt-L5-P1","Level 5","");

test_table ("adapt-L5-P1-mesh",
	    array("0.000000","0.020000","0.040000","0.060000",
		  "0.080000","0.100000"), $types);
test_table ("adapt-L5-P1-de",
	    array("0.000000","0.020000","0.040000","0.060000",
		  "0.080000","0.100000"), $types);
test_table ("adapt-L5-P1-te",
	    array("0.000000","0.020000","0.040000","0.060000",
		  "0.080000","0.100000"), $types);
test_table ("adapt-L5-P1-vx",
	    array("0.000000","0.020000","0.040000","0.060000",
		  "0.080000","0.100000"), $types);
test_table ("adapt-L5-P1-vy",
	    array("0.000000","0.020000","0.040000","0.060000",
		  "0.080000","0.100000"), $types);

end_hidden("adapt_L5");


//======================================================================

test_group("Balance");

begin_hidden("balance_none", "None");

tests("Enzo","enzo-p","test_balance_none","None","Balance");
test_table ("Balance/None/balance-mesh",
	    array("00000","00002","00004","00006","00008","00010","00020"), $types);
test_table ("Balance/None/balance-de",
	    array("00000","00002","00004","00006","00008","00010","00020"), $types);

end_hidden("balance_none");

begin_hidden("balance_rotate", "RotateLB");

tests("Enzo","enzo-p","test_balance_rotate","Rotate","Balance");
test_table ("Balance/Rotate/balance-mesh",
	    array("00000","00002","00004","00006","00008","00010","00020"), $types);
test_table ("Balance/Rotate/balance-de",
	    array("00000","00002","00004","00006","00008","00010","00020"), $types);

end_hidden("balance_rotate");

begin_hidden("balance_greedy", "GreedyLB");

tests("Enzo","enzo-p","test_balance_greedy","Greedy","Balance");
test_table ("Balance/Greedy/balance-mesh",
	    array("00000","00002","00004","00006","00008","00010","00020"), $types);
test_table ("Balance/Greedy/balance-de",
	    array("00000","00002","00004","00006","00008","00010","00020"), $types);

end_hidden("balance_greedy");

begin_hidden("balance_hybrid", "HybridLB");

tests("Enzo","enzo-p","test_balance_hybrid","Hybrid","Balance");
test_table ("Balance/Hybrid/balance-mesh",
	    array("00000","00002","00004","00006","00008","00010","00020"), $types);
test_table ("Balance/Hybrid/balance-de",
	    array("00000","00002","00004","00006","00008","00010","00020"), $types);

end_hidden("balance_hybrid");

// begin_hidden("balance_neighbor", "NeighborLB");

// tests("Enzo","enzo-p","test_balance_neighbor","Neighbor","Balance");
// test_table ("Balance/Neighbor/balance-mesh",
//	    array("00000","00002","00004","00006","00008","00010","00020"), $types);
//test_table ("Balance/Neighbor/balance-de",
//	    array("00000","00002","00004","00006","00008","00010","00020"), $types);

//end_hidden("balance_neighbor");

begin_hidden("balance_rand_cent", "RandCentLB");

tests("Enzo","enzo-p","test_balance_rand_cent","RandCent","Balance");
test_table ("Balance/RandCent/balance-mesh",
	    array("00000","00002","00004","00006","00008","00010","00020"), $types);
test_table ("Balance/RandCent/balance-de",
	    array("00000","00002","00004","00006","00008","00010","00020"), $types);

end_hidden("balance_rand_cent");

begin_hidden("balance_refine", "RefineLB");

tests("Enzo","enzo-p","test_balance_refine","Refine","Balance");
test_table ("Balance/Refine/balance-mesh",
	    array("00000","00002","00004","00006","00008","00010","00020"), $types);
test_table ("Balance/Refine/balance-de",
	    array("00000","00002","00004","00006","00008","00010","00020"), $types);

end_hidden("balance_refine");

//======================================================================

test_group("Boundary");

begin_hidden("boundary_reflecting_2d", "2D Reflecting");

tests("Enzo","enzo-p","test_boundary_reflecting-2d","Reflecting 2D","");
test_table ("boundary_reflecting-2d",
	    array("0000","0100","0200","0300","0400"), $types);
end_hidden("boundary_reflecting_2d");

begin_hidden("boundary_periodic_2d", "2D Periodic");

tests("Enzo","enzo-p","test_boundary_periodic-2d","Periodic 2D","");
test_table ("boundary_periodic-2d",
	    array("0000","0100","0200","0300","0400"), $types);
end_hidden("boundary_periodic_2d");

begin_hidden("boundary_outflow_2d", "2D Outflow");

tests("Enzo","enzo-p","test_boundary_outflow-2d","Outflow 2D","");
test_table ("boundary_outflow-2d",
	    array("0000","0100","0200","0300","0400"), $types);
end_hidden("boundary_outflow_2d");

//----------------------------------------------------------------------

begin_hidden("boundary_reflecting_3d", "3D Reflecting");

tests("Enzo","enzo-p","test_boundary_reflecting-3d","Reflecting 3D","");
test_table ("boundary_reflecting-3d",
	    array("0000","0020","0040","0060","0080"), $types);
end_hidden("boundary_reflecting_3d");

begin_hidden("boundary_periodic_3d", "3D Periodic");

tests("Enzo", "enzo-p","test_boundary_periodic-3d","Periodic 3D","");
test_table ("boundary_periodic-3d",
	    array("0000","0020","0040","0060","0080"), $types);
end_hidden("boundary_periodic_3d");

begin_hidden("boundary_outflow_3d", "3D Outflow");

tests("Enzo","enzo-p","test_boundary_outflow-3d","Outflow 3D","");
test_table ("boundary_outflow-3d",
	    array("0000","0020","0040","0060","0080"), $types);
end_hidden("boundary_outflow_3d");

//----------------------------------------------------------------------

test_group("Initial");

begin_hidden("initial_mask","png mask initial conditions");

tests("Enzo","enzo-p","test_initial_png","","");
test_table ("initial_png",
	    array("00","10","20","30","40", "50"), $types);
end_hidden("initial_mask");

//----------------------------------------------------------------------

test_group("Output");

begin_hidden("output_stride_1", "Stride 1");
tests("Enzo","enzo-p","test_output-stride-1","","");
test_table_blocks ("output-stride-1", array("00","10","20"),$types);
end_hidden("output_stride_1");

begin_hidden("output_stride_2", "Stride 2");
tests("Enzo","enzo-p","test_output-stride-2","","");
test_table_blocks ("output-stride-2",  array("00","10","20"), $types);
end_hidden("output_stride_2");

begin_hidden("output_stride_3", "Stride 3");
tests("Enzo","enzo-p","test_output-stride-4","","");
test_table_blocks ("output-stride-4",  array("00","10","20"), $types);
end_hidden("output_stride_3");

//======================================================================

     echo "<a name=\"cello\"><h1>Cello unit tests</h1>";

test_group("Disk");

begin_hidden("disk_hdf5", "HDF5");
tests("Cello","test_FileHdf5", "test_FileHdf5","","");
end_hidden("disk_hdf5");
begin_hidden("disk_ifrit", "IFRIT");
tests("Cello","test_FileIfrit","test_FileIfrit","","");
end_hidden("disk_ifrit");

//----------------------------------------------------------------------

test_group("Error");

begin_hidden("error", "Error");
tests("Cello","test_Error","test_Error","","");
end_hidden("error");


//----------------------------------------------------------------------

test_group("Field");

begin_hidden("field", "Field");
tests("Cello","test_Field","test_Field","","");
end_hidden("field");
begin_hidden("field_descr", "FieldDescr");
tests("Cello","test_FieldDescr","test_FieldDescr","","");
end_hidden("field_descr");
begin_hidden("field_data", "FieldData");
tests("Cello","test_FieldData","test_FieldData","","");
end_hidden("field_data");
begin_hidden("field_face", "FieldFace");
tests("Cello","test_FieldFace","test_FieldFace","","");
end_hidden("field_face");
begin_hidden("it_field", "ItField");
tests("Cello","test_ItField","test_ItField","","");
end_hidden("it_field");
begin_hidden("grouping", "Grouping");
tests("Cello","test_Grouping","test_Grouping","","");
end_hidden("grouping");

//----------------------------------------------------------------------

test_group("Memory");


begin_hidden("memory", "Memory");
tests("Cello","test_Memory","test_Memory","","");
end_hidden("memory");


//----------------------------------------------------------------------

test_group("Mesh");

begin_hidden("hierarchy", "Hierarchy");
tests("Cello","test_Hierarchy","test_Hierarchy","",""); 
end_hidden("hierarchy");
begin_hidden("data", "Data");
tests("Cello","test_Data","test_Data","",""); 
end_hidden("data");
begin_hidden("index", "Index");
tests("Cello","test_Index","test_Index","",""); 
end_hidden("index");
begin_hidden("tree", "Tree");
tests("Cello","test_Tree","test_Tree","",""); 
end_hidden("tree");
begin_hidden("it_face", "ItFace");
tests("Cello","test_ItFace","test_ItFace","",""); 
end_hidden("it_face");

begin_hidden("tree_initial", "Tree (initial)");
printf ("<img width=257 src=\"test_tree_1-initial.png\"></img></br>\n");
end_hidden("tree_initial");
begin_hidden("tree_balanced", "Tree (balanced)");
printf ("<img width=257 src=\"test_tree_2-balanced.png\"></img>\n");
end_hidden("tree_balanced");
begin_hidden("tree_merged", "Tree (merged)");
printf ("<img width=257 src=\"test_tree_3-merged.png\"></img></br>\n");
end_hidden("tree_merged");

begin_hidden("tree_density", "TreeDensity");
tests("Cello","test_TreeDensity","test_TreeDensity","",""); 
end_hidden("tree_density");

begin_hidden("mesh_density_2d", "Density 2D");
printf ("<img width=257 src=\"density_xy_1-initial.png\"></img>\n");
printf ("<img width=257 src=\"density_xy_2-balanced.png\"></img>\n");
printf ("<img width=257 src=\"density_xy_3-coalesced.png\"></img></br>\n");
end_hidden("mesh_density_coalesced_2d");

begin_hidden("mesh_density_3d", "Density 3D");
printf ("<img width=257 src=\"density_3d_1-initial.png\"></img>\n");
printf ("<img width=257 src=\"density_3d_2-balanced.png\"></img>\n");
printf ("<img width=257 src=\"density_3d_3-coalesced.png\"></img></br>\n");
end_hidden("mesh_density_coalesced_3d");

begin_hidden("node", "Node");
tests("Cello","test_Node","test_Node","",""); 
end_hidden("node");
begin_hidden("node_trace", "NodeTrace");
tests("Cello","test_NodeTrace","test_NodeTrace","",""); 
end_hidden("node_trace");
begin_hidden("it_node", "ItNode");
tests("Cello","test_ItNode","test_ItNode","",""); 
end_hidden("it_node");

//----------------------------------------------------------------------

test_group("Monitor");

begin_hidden("monitor", "Monitor");
tests("Cello","test_Monitor","test_Monitor","","");
end_hidden("monitor");

//----------------------------------------------------------------------

test_group("Parameters");

begin_hidden("parameters", "Parameters");
tests("Cello","test_Parameters","test_Parameters","","");
end_hidden("parameters");

//----------------------------------------------------------------------

test_group("Particle");

begin_hidden("particle", "Particle");
tests("Cello","test_Particle","test_Particle","","");
end_hidden("particle");
//----------------------------------------------------------------------

test_group("Performance");

begin_hidden("performance", "Performance");
tests("Cello","test_Performance","test_Performance","","");
end_hidden("performance");
begin_hidden("papi", "Papi");
tests("Cello","test_Papi",       "test_Papi","","");
end_hidden("papi");
begin_hidden("timer", "Timer");
tests("Cello","test_Timer",       "test_Timer","","");
end_hidden("timer");

//----------------------------------------------------------------------

test_group("Problem");

begin_hidden("mask", "Mask");
tests("Cello","test_Mask",   "test_Mask","","");
end_hidden("mask");
begin_hidden("refresh", "Refresh");
tests("Cello","test_Refresh","test_Refresh","","");
end_hidden("refresh");
begin_hidden("value", "Value");
tests("Cello","test_Value",   "test_Value","","");
end_hidden("value");

//----------------------------------------------------------------------

test_group("Schedule");

begin_hidden("schedule", "Schedule");
tests("Cello","test_Schedule","test_Schedule","","");
end_hidden("schedule");

//----------------------------------------------------------------------

test_group("Colormap");

begin_hidden("colormap", "Colormap");
tests("Cello","test_Colormap","test_Colormap","","");
end_hidden("colormap");

?>
</br/>
</body>
</html>
