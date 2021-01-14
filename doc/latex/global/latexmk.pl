
## Defaults
$do_cd    = 1;   ## Change to the directory of the main source file
$silent   = 1;   ## Less verbosity

## Use of 'build' relative directory
$ENV{'openout_any'}='a';
$out_dir = '../build';

## Global option
set_tex_cmds( '-shell-escape' );

$makeindex = "makeindex %O -s ../../global/index -o $out_dir/%D $out_dir/%S";
