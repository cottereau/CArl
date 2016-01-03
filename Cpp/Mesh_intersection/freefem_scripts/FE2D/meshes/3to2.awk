## Author: G. D. McBain <gmcbain>
## Created: 2009-09-08
## Keywords: meshing, two-dimensional, MEDIT .mesh

## Usage: gawk -f 3to2.awk

## Input: A MEDIT .mesh file, as for example output by Gmsh.

## Output: A modified MEDIT .mesh file with the "Dimension" field set
## to 2 and the z-coordinate stripped from each vertex.  This is the
## format FreeFem++ wants for its readmesh command.

## Examples: msmdir.000708

BEGIN { 
  dimension = 0;		# flag for processing Dimension
  vertices = 0;			# number of vertices to process
}
dimension == 0 && vertices <= 0 { print } # normal (pass unmodified)

dimension == -1 { 		# Output Dimension as 2.
  print 2;
  dimension = 0;		# our work is done
}
/Dimension/ { dimension = -1; }	# Get ready to process Dimension.

vertices > 0 {			# Omit z for each vertex.
  print $1, $2, $4;
  vertices--;
}
vertices == -1 { vertices = $1; } # Read the number of vertices.
/Vertices/ { vertices = -1; }	# Get ready to process vertices.
