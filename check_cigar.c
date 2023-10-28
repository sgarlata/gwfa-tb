#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
	FILE *fcigar;
	char edit;
	int count;
	int match = 0, mismatch = 0, deletion = 0, insertion = 0;

	fcigar = fopen(argv[1], "r");

	if (fcigar == NULL)
	{
		fprintf(stderr, "Error opening file.\n");
		return 1;
	}

	while (fscanf(fcigar, "%d%c", &count, &edit) != EOF)
	{
		switch (edit)
		{
		case '=':
			match += count;
			break;

		case 'X':
			mismatch += count;
			break;

		case 'D':
			deletion += count;
			break;

		case 'I':
			insertion += count;
			break;

		default:
			fprintf(stderr, "Edits error.\n");
			break;
		}
	}
	fprintf(stdout, "Edit distance: %d\n", mismatch + deletion + insertion);
	fprintf(stdout, "N. of matches: %d\n", match);
	fprintf(stdout, "N. of mismatches: %d\n", mismatch);
	fprintf(stdout, "N. of deletions: %d\n", deletion);
	fprintf(stdout, "N. of insertions: %d\n", insertion);

	fclose(fcigar);
	return 0;
}
