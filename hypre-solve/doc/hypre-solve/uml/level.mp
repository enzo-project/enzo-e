input metauml;

beginfig(1);

   ClassName.level("Level");
   ClassName.levelstructured("LevelStructured");
   ClassName.levelunstructured("LevelUnstructured");

   centered_align_top(levelstructured, levelunstructured)(30, below(level.s, 30));

   drawObjects(level,levelstructured,levelunstructured);

   link(inheritance)(pathStepY(levelstructured.n, level.s, 10));
   link(inheritance)(pathStepY(levelunstructured.n, level.s, 10));

endfig;
end

   