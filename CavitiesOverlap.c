#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdint.h>

// Windows or Linux specific functions for bit count
#if defined(_WIN64)
#include <intrin.h>
#define CountBitsSet (unsigned int)__popcnt64
#elif defined(__x86_64__)
#define CountBitsSet (unsigned int)__builtin_popcountll
#endif

// Windows or Linux specific libraries for filesystem/directories processing
#if defined(_WIN32) || defined(_WIN64)
	#define WIN32_LEAN_AND_MEAN
	#include <Windows.h>
	#define DIR_SEP '\\'
#elif defined(__unix__) || defined(__linux__)
	#include <dirent.h>
	#define DIR_SEP '/'
#endif

#define PAUSE_AT_EXIT 1
#define ARRAY_GROW_SIZE 256
#define MAX_SPHERE_RADIUS 1.0f
#define GRID_SIZE 0.1f

#define FILENAME_STATS "_ALL-stats.csv"
#define FILENAME_ERRORS "_ALL-errors.csv"

typedef struct _SphereCoords {
	float x;
	float y;
	float z;
	//float r;
	//float intersected;
} SphereCoords;

typedef struct _CavityInfo {
	int numSpheres;
	int maxNumSpheres;
	SphereCoords *spheres;
	double totalVolume;
	double totalOverlap;
} CavityInfo;

int numGTCavities = 0;
int numMRCavities = 0;
CavityInfo *gtCavities = NULL;
CavityInfo *mrCavities = NULL;

FILE *statsFile = NULL;
FILE *errorsFile = NULL;
char *outputDir = NULL;

// Create new file/path name string
// Option 'mode':
//  0 - full path name ,
//  1 - remove file extension ,
//  2 - keep only file name (remove path) ,
//  3 - keep only path (remove file name; add extra slash if needed) ,
//  4 - keep only extension (remove path and base name; keep the dot)
//  5 - keep only base name (remove path and file extension)
// Option 'changecase':
//  0 - do not change case
//  1 - convert to lowercase
//  2 - convert to uppercase
char* CreateFilename(char *filename, char *extrastart, char *extraend, int mode, int changecase) {
	char *charptr, *resultfilename;
	int startpos, endpos, len, pos;
	startpos = (-1);
	endpos = 0;
	len = 0;
	charptr = filename;
	while ((*charptr) != '\0') {
		if ((*charptr) == DIR_SEP) startpos = len; // set start to position of last slash character
		if ((*charptr) == '.') endpos = len; // set end to position of last dot character
		charptr++;
		len++; // string length of filename
	}
	if (endpos == 0 || endpos < startpos){ // if no dot character was found or if the dot appears before the slash, the full string is a directory, not a file
		endpos = len; // simulate a non-existent dot at the end
		if (startpos != (len - 1)) startpos = len; // if the slash is not already the last char, set a ficticious slash at the end of the string
	} // if a slash was not found, it is already correctly set at (-1) if the string is a file name (with a dot), or set above at (len) if the string is a directory (without a dot)
	if (mode == 0) { // do not cut anything ; set to beginning and end of full string
		startpos = (-1);
		endpos = len;
	}
	else if (mode == 1) { // cut extension only
		startpos = (-1);
	}
	else if (mode == 2) { // cut path only (keep only file name)
		endpos = len;
	}
	else if (mode == 3) { // cut file name (keep only path)
		if (startpos == (-1)) endpos = 0; // no slash was found but a dot was (otherwise 'startpos' would have been changed above to 'len'), so the full string is a file name
		else endpos = startpos;
		startpos = (-1);
	}
	else if (mode == 4) { // cut both path and base name (keep only extension)
		startpos = endpos;
		startpos--; // keep the dot character
		endpos = len;
	} // else cut both path and extension (mode==5) (keep only base name)
	len = 0; // string length of extra start and end strings
	if (extrastart != NULL) {
		for (charptr = extrastart; (*charptr) != '\0'; len++, charptr++);
	}
	if (extraend != NULL) {
		for (charptr = extraend; (*charptr) != '\0'; len++, charptr++);
	}
	if (mode == 3) {
		if (endpos == 0) { // if we want the path but the resulting string was empty, add an initial "./"
			len++;
			if (extrastart != NULL) len++;
		}
		len++; // add an extra slash at the end
	}
	resultfilename = (char *)calloc(((endpos - startpos - 1) + len + 1), sizeof(char));
	len = 0; // current position in the new string
	if (mode == 3 && endpos == 0) { // add an initial "./" to the path if needed
		resultfilename[len++] = '.';
		if (extrastart != NULL) resultfilename[len++] = DIR_SEP;
	}
	if (extrastart != NULL) { // copy extra start chars
		for (charptr = extrastart; (*charptr) != '\0'; resultfilename[len] = (*charptr), len++, charptr++);
	}
	if (changecase == 0) { // copy middle chars 
		for (pos = (startpos + 1); pos < endpos; resultfilename[len++] = filename[pos++]);
	}
	else { // copy middle chars but convert them to lowercase/uppercase
		charptr = (filename + startpos + 1);
		for (pos = (startpos + 1); pos < endpos; charptr++, pos++, len++) {
			resultfilename[len] = (*charptr);
			if ((changecase == 1) && (*charptr) >= 'A' && (*charptr) <= 'Z') resultfilename[len] += 32;
			else if ((changecase == 2) && (*charptr) >= 'a' && (*charptr) <= 'z') resultfilename[len] -= 32;
		}
	}
	if (mode == 3) { // add a final dash to the path
		resultfilename[len++] = DIR_SEP;
	}
	if (extraend != NULL) { // copy extra end chars
		for (charptr = extraend; (*charptr) != '\0'; resultfilename[len] = (*charptr), len++, charptr++);
	}
	resultfilename[len] = '\0';
	return resultfilename;
}

void FreeCavityArrays() {
	int i;
	if (gtCavities != NULL) {
		for (i = 0; i < numGTCavities; i++) free(gtCavities[i].spheres);
		free(gtCavities);
		gtCavities = NULL;
		numGTCavities = 0;
	}
	if (mrCavities != NULL) {
		for (i = 0; i < numMRCavities; i++) free(mrCavities[i].spheres);
		free(mrCavities);
		mrCavities = NULL;
		numMRCavities = 0;
	}
}

// Returns the number of cavities loaded or the line number (negative) where an error occured
int LoadCavitiesFile(char *cavFilename, CavityInfo **cavities) {
	FILE *cavFile;
	int lineType, numLines, cavityId, maxCavityId, numValuesRead, i;
	//int prevCavityId;
	unsigned int totalNumSpheres;
	float x, y, z;
	char c, nl[2];
	char *tempFilename;
#ifdef _WIN32
	while (1) {
		tempFilename = cavFilename;
		if ((cavFile = fopen(tempFilename, "r")) != NULL) break;
		tempFilename = CreateFilename(cavFilename, NULL, ".csv", 1, 0); // try with .csv extension
		if ((cavFile = fopen(tempFilename, "r")) != NULL) break;
		free(tempFilename);
		printf("\n> ERROR: Cavities file \"%s\" not found\n", cavFilename);
		return 0;
	}
#else
	int n;
	while (1) {
		tempFilename = CreateFilename(cavFilename, NULL, NULL, 0, 0);
		for (n = 0; tempFilename[n] != '\0'; n++); // go to end of string
		for (i = 0; i < 2; i++) { // try two extensions (.txt and .csv) in both lettercases
			if ((cavFile = fopen(tempFilename, "r")) != NULL) break; // try opening regular case file name
			while (n >= 0 && tempFilename[n] != DIR_SEP) n--; // go to position of last dash
			for (n++; (c = tempFilename[n]) != '.' && c != '\0'; n++) if (c >= 'A' && c <= 'Z') tempFilename[n] = (c + 32); // change case of base name before the extension
			if ((cavFile = fopen(tempFilename, "r")) != NULL) break; // try opening lowercase file name
			while (n >= 0 && tempFilename[n] != DIR_SEP) n--;
			for (n++; (c = tempFilename[n]) != '.' && c != '\0'; n++) if (c >= 'a' && c <= 'z') tempFilename[n] = (c - 32);
			if ((cavFile = fopen(tempFilename, "r")) != NULL) break; // try opening uppercase file name
			if (i == 0) { // change file extension
				while (n >= 0 && tempFilename[n] != DIR_SEP) n--;
				for (n++; (c = tempFilename[n]) != '.' && c != '\0'; n++) tempFilename[n] = cavFilename[n]; // back to original case
				if (tempFilename[++n] != '\0') tempFilename[n] = 'c';
				if (tempFilename[++n] != '\0') tempFilename[n] = 's';
				if (tempFilename[++n] != '\0') tempFilename[n] = 'v';
			}
		}
		if (i != 2) break; // file was successfully opened
		free(tempFilename);
		printf("\n> ERROR: Cavities file \"%s\" not found\n", cavFilename);
		return 0;
	}
#endif
	if (tempFilename != cavFilename) { // save new file name if it was changed
		for (i = 0; tempFilename[i] != '\0'; i++) cavFilename[i] = tempFilename[i];
		free(tempFilename);
	}
	//printf("<%s> ", cavFilename);
	lineType = (-1); // check the correct line format
	i = 0; // count of dot characters
	while ((c = fgetc(cavFile)) != '\n' && c != EOF) {
		if (c == '.') i++;
		if (c == ' ' || c == '\t') {
			if (i == 0) lineType = 0; // first number does not have a dot (integer at 1st position)
			else if (i == 3) lineType = 1; // three numbers with dots seen before (integer at 4th position)
		}
	}
	if (c == EOF) {
		printf("\n> ERROR: No data found inside file \"%s\"\n", cavFilename);
		return 0;
	}
	if (i != 3 || lineType == (-1)) {
		printf("\n> ERROR: Invalid line format in file \"%s\"\n", cavFilename);
		return (-1);
	}
	rewind(cavFile);
	numLines = 0;
	totalNumSpheres = 0;
	cavityId = (-1);
	//prevCavityId = (-1);
	maxCavityId = (-1);
	while (1) { // load cavities and corresponding spheres
		numLines++;
		if (lineType == 0) numValuesRead = fscanf(cavFile, "%d %f %f %f%1[\n]", &cavityId, &x, &y, &z, nl); // read extra '\n' char explicitly at end to check for line format consistency
		else numValuesRead = fscanf(cavFile, "%f %f %f %d%1[\n]", &x, &y, &z, &cavityId, nl);
		if (numValuesRead != 5) break; // check end of file
		if (cavityId<0 || cavityId>(maxCavityId + 100)) break; // invalid cavity id
		if (cavityId > maxCavityId) { // new cavity
			(*cavities) = (CavityInfo *)realloc((*cavities), (cavityId + 1) * sizeof(CavityInfo));
			for (i = (maxCavityId + 1); i <= cavityId; i++) {
				(*cavities)[i].spheres = (SphereCoords *)calloc(ARRAY_GROW_SIZE, sizeof(SphereCoords));
				(*cavities)[i].maxNumSpheres = ARRAY_GROW_SIZE;
				(*cavities)[i].numSpheres = 0;
				(*cavities)[i].totalVolume = 0.0;
				(*cavities)[i].totalOverlap = 0.0;
			}
			maxCavityId = cavityId;
		}
		if (((*cavities)[cavityId].numSpheres) == ((*cavities)[cavityId].maxNumSpheres)) {
			(*cavities)[cavityId].maxNumSpheres += ARRAY_GROW_SIZE;
			(*cavities)[cavityId].spheres = (SphereCoords *)realloc((*cavities)[cavityId].spheres, ((*cavities)[cavityId].maxNumSpheres) * sizeof(SphereCoords));
		}
		i = (*cavities)[cavityId].numSpheres;
		(*cavities)[cavityId].spheres[i].x = x;
		(*cavities)[cavityId].spheres[i].y = y;
		(*cavities)[cavityId].spheres[i].z = z;
		//(*cavities)[cavityId].spheres[i].intersected = 0.0;
		(*cavities)[cavityId].numSpheres++;
		totalNumSpheres++;
		//prevCavityId = cavityId;
	}
	fclose(cavFile);
	if (numValuesRead != EOF) { // line with invalid format
		if (cavityId<0 || cavityId>(maxCavityId + 100)) printf("\n> ERROR: Invalid cavity id %d in line %d of file \"%s\"\n", cavityId, numLines, cavFilename);
		else printf("\n> ERROR: Invalid format in line %d of file \"%s\"\n", numLines, cavFilename);
		return (-numLines);
	}
	if (maxCavityId == (-1)) {
		printf("\n> ERROR: No cavities found in file \"%s\"\n", cavFilename);
		return (-numLines);
	}
	maxCavityId++; // number of cavities
	printf("(%d)(%u) ", maxCavityId, totalNumSpheres);
	return maxCavityId;
}

void LoadCavityIntoGrid(CavityInfo *cavity, uint64_t *grid, float gridSize, unsigned int gridSizeXYZ[3], float minXYZ[3]) {
	unsigned int numCavitySpheres, sphereId, sphereRadius;
	unsigned int cx, cy, cz, x, y, z, dx, dy, dz;
	unsigned int gridPoint;
	SphereCoords *sphere;
	sphereRadius = (unsigned int)floorf((float)MAX_SPHERE_RADIUS / gridSize); // sphere radius in grid points
	numCavitySpheres = (cavity->numSpheres);
	for (sphereId = 0; sphereId < numCavitySpheres; sphereId++) { // process all spheres of this cavity
		sphere = &(cavity->spheres[sphereId]);
		cx = (unsigned int)floorf(((sphere->x) - minXYZ[0]) / gridSize); // convert position to grid coordinates
		cy = (unsigned int)floorf(((sphere->y) - minXYZ[1]) / gridSize);
		cz = (unsigned int)floorf(((sphere->z) - minXYZ[2]) / gridSize);
		dx = sphereRadius;
		for (x = 0; x <= dx; x++) { // set all points of this sphere
			dy = (unsigned int)floorf(sqrtf((float)(sphereRadius*sphereRadius - x * x))); // y = V(R^2-x^2)
			for (y = 0; y <= dy; y++) {
				dz = (unsigned int)floorf(sqrtf((float)(sphereRadius*sphereRadius - x * x - y * y))); // z = V(R^2-x^2-y^2)
				for (z = 0; z <= dz; z++) { // grid point = z*(SX*SY) + y*SX + x ; set all 8 points at once, one for each octant
					gridPoint = ((cz + z)*gridSizeXYZ[1] + (cy + y))*gridSizeXYZ[0] + (cx + x); // (+X,+Y,+Z)
					grid[(gridPoint >> 6)] |= (1ULL << (gridPoint & 0x3F)); // set bit to 1 in grid ; wordId=(point/64) ; wordOffset=(point%64)
					gridPoint = ((cz + z)*gridSizeXYZ[1] + (cy + y))*gridSizeXYZ[0] + (cx - x); // (-X,+Y,+Z)
					grid[(gridPoint >> 6)] |= (1ULL << (gridPoint & 0x3F));
					gridPoint = ((cz + z)*gridSizeXYZ[1] + (cy - y))*gridSizeXYZ[0] + (cx + x); // (+X,-Y,+Z)
					grid[(gridPoint >> 6)] |= (1ULL << (gridPoint & 0x3F));
					gridPoint = ((cz + z)*gridSizeXYZ[1] + (cy - y))*gridSizeXYZ[0] + (cx - x); // (-X,-Y,+Z)
					grid[(gridPoint >> 6)] |= (1ULL << (gridPoint & 0x3F));
					gridPoint = ((cz - z)*gridSizeXYZ[1] + (cy + y))*gridSizeXYZ[0] + (cx + x); // (+X,+Y,-Z)
					grid[(gridPoint >> 6)] |= (1ULL << (gridPoint & 0x3F));
					gridPoint = ((cz - z)*gridSizeXYZ[1] + (cy + y))*gridSizeXYZ[0] + (cx - x); // (-X,+Y,-Z)
					grid[(gridPoint >> 6)] |= (1ULL << (gridPoint & 0x3F));
					gridPoint = ((cz - z)*gridSizeXYZ[1] + (cy - y))*gridSizeXYZ[0] + (cx + x); // (+X,-Y,-Z)
					grid[(gridPoint >> 6)] |= (1ULL << (gridPoint & 0x3F));
					gridPoint = ((cz - z)*gridSizeXYZ[1] + (cy - y))*gridSizeXYZ[0] + (cx - x); // (-X,-Y,-Z)
					grid[(gridPoint >> 6)] |= (1ULL << (gridPoint & 0x3F));
				}
			}
		}
	} // loop for all spheres of this cavity
}

#if !defined(_WIN64) && !defined(__x86_64__)
inline unsigned int CountBitsSet(uint64_t word) {
	unsigned int count;
	for (count = 0; word != 0ULL; count++) word &= (word - 1ULL); // keep clearing the rightmost set bit while there are still bits set
	return count;
}
#endif

int ProcessCavitiesWithGrid(float gridSize, char *filename, int createVis) {
	char *outputFilename, *visualOutputFilename, *fileId;
	FILE *outputFile, *visualOutputFile;
	int gtCavId, mrCavId;
	unsigned int numCavitySpheres, sphereId;
	SphereCoords *sphere;
	int numTP, numFP, numFN;
	double pctTP, pctFP, pctFN, pctTN, sumOverlaps, gridPointVolume;
	unsigned int totalMRCavitiesVolume, totalGTCavitiesVolume;
	float minX, minY, minZ, minCoords[3];
	float maxX, maxY, maxZ;
	float rangeX, rangeY, rangeZ;
	unsigned int gridSizeX, gridSizeY, gridSizeZ, gridSizes[3];
	unsigned int numBitWords, wordId, numSetGridPoints, gridPoint;
	uint64_t *gtBitWordsGrid, *mrBitWordsGrid, tempBitWord;
	float visualGridSize, x, y, z;
	ldiv_t divres;
	minX = FLT_MAX;
	minY = FLT_MAX;
	minZ = FLT_MAX;
	maxX = -(FLT_MAX);
	maxY = -(FLT_MAX);
	maxZ = -(FLT_MAX);
	for (gtCavId = 0; gtCavId < numGTCavities; gtCavId++) { // loop through all spheres of all GT cavities
		numCavitySpheres = gtCavities[gtCavId].numSpheres;
		for (sphereId = 0; sphereId < numCavitySpheres; sphereId++) {
			sphere = &(gtCavities[gtCavId].spheres[sphereId]);
			if ((sphere->x) < minX) minX = (sphere->x);
			if ((sphere->x) > maxX) maxX = (sphere->x);
			if ((sphere->y) < minY) minY = (sphere->y);
			if ((sphere->y) > maxY) maxY = (sphere->y);
			if ((sphere->z) < minZ) minZ = (sphere->z);
			if ((sphere->z) > maxZ) maxZ = (sphere->z);
		}
	}
	for (mrCavId = 0; mrCavId < numMRCavities; mrCavId++) { // loop through all spheres of all MR cavities
		numCavitySpheres = mrCavities[mrCavId].numSpheres;
		for (sphereId = 0; sphereId < numCavitySpheres; sphereId++) {
			sphere = &(mrCavities[mrCavId].spheres[sphereId]);
			if ((sphere->x) < minX) minX = (sphere->x);
			if ((sphere->x) > maxX) maxX = (sphere->x);
			if ((sphere->y) < minY) minY = (sphere->y);
			if ((sphere->y) > maxY) maxY = (sphere->y);
			if ((sphere->z) < minZ) minZ = (sphere->z);
			if ((sphere->z) > maxZ) maxZ = (sphere->z);
		}
	}
	minX -= (MAX_SPHERE_RADIUS + gridSize);
	maxX += (MAX_SPHERE_RADIUS + gridSize);
	minY -= (MAX_SPHERE_RADIUS + gridSize);
	maxY += (MAX_SPHERE_RADIUS + gridSize);
	minZ -= (MAX_SPHERE_RADIUS + gridSize);
	maxZ += (MAX_SPHERE_RADIUS + gridSize);
	rangeX = ceilf(maxX - minX);
	rangeY = ceilf(maxY - minY);
	rangeZ = ceilf(maxZ - minZ);
	gridSizeX = (unsigned int)ceilf(rangeX / gridSize); // number of grid points (bits) ; +1 to account for floor/ceil roundings
	gridSizeY = (unsigned int)ceilf(rangeY / gridSize);
	gridSizeZ = (unsigned int)ceilf(rangeZ / gridSize);
	minCoords[0] = minX; // arrays of variables for LoadCavityIntoGrid function
	minCoords[1] = minY;
	minCoords[2] = minZ;
	gridSizes[0] = gridSizeX;
	gridSizes[1] = gridSizeY;
	gridSizes[2] = gridSizeZ;
	printf("(%ux%ux%u)", gridSizeX, gridSizeY, gridSizeZ);
	numBitWords = ((gridSizeX * gridSizeY * gridSizeZ) + 63) / 64; // total number of 64bit words required
	if (numBitWords < (gridSizeX + gridSizeY + gridSizeZ)) { // check 32bit word overflow
		printf("\n> ERROR: Grid size is too large\n");
		return (-1);
	}
	printf("(%uMB) ", (unsigned int)(numBitWords * sizeof(uint64_t) / (1024 * 1024)));
	gtBitWordsGrid = (uint64_t *)calloc(numBitWords, sizeof(uint64_t));
	mrBitWordsGrid = (uint64_t *)calloc(numBitWords, sizeof(uint64_t));
	if (gtBitWordsGrid == NULL || mrBitWordsGrid == NULL) {
		printf("\n> ERROR: Not enough memory to create grids\n");
		return (-1);
	}
	outputFilename = CreateFilename(filename, outputDir, "-output.txt", 5, 1);
	if ((outputFile = fopen(outputFilename, "w")) == NULL) {
		printf("\n> ERROR: Cannot create output file \"%s\"\n", outputFilename);
		return (-1);
	}
	free(outputFilename);
	fileId = CreateFilename(filename, NULL, NULL, 5, 2);
	fprintf(outputFile, "%s", fileId);
	for (gtCavId = 0; gtCavId < numGTCavities; gtCavId++) fprintf(outputFile, " %d", gtCavId);
	fprintf(outputFile, "\n");
	totalMRCavitiesVolume = 0;
	totalGTCavitiesVolume = 0;
	gridPointVolume = (double)pow((double)gridSize, 3.0); // volume of each grid voxel
	for (mrCavId = 0; mrCavId < numMRCavities; mrCavId++) { // process all MR cavities
		LoadCavityIntoGrid(&(mrCavities[mrCavId]), mrBitWordsGrid, gridSize, gridSizes, minCoords); // load MR cavity into MR grid
		numSetGridPoints = 0;
		for (wordId = 0; wordId < numBitWords; wordId++) {
			if (mrBitWordsGrid[wordId] == (uint64_t)0) continue;
			numSetGridPoints += CountBitsSet(mrBitWordsGrid[wordId]); // count number of bits set to 1 in a 64-bit word
		}
		totalMRCavitiesVolume += numSetGridPoints;
		mrCavities[mrCavId].totalVolume = (double)numSetGridPoints;
		mrCavities[mrCavId].totalOverlap = 0.0;
		fprintf(outputFile, "%d", mrCavId);
		for (gtCavId = 0; gtCavId < numGTCavities; gtCavId++) { // for each MR cavity, compute overlap with each GT cavity
			LoadCavityIntoGrid(&(gtCavities[gtCavId]), gtBitWordsGrid, gridSize, gridSizes, minCoords); // load GT cavity into GT grid
			if (mrCavId == 0) { // get total volume of this GT cavity, but only one time, when processing 1st MR cavity
				numSetGridPoints = 0;
				for (wordId = 0; wordId < numBitWords; wordId++) {
					if (gtBitWordsGrid[wordId] == (uint64_t)0) continue;
					numSetGridPoints += CountBitsSet(gtBitWordsGrid[wordId]); // count number of bits set to 1 in a 64-bit word
				}
				totalGTCavitiesVolume += numSetGridPoints;
				//gtCavities[gtCavId].totalVolume = ((double)numSetGridPoints * gridPointVolume);
				gtCavities[gtCavId].totalVolume = (double)numSetGridPoints;
				gtCavities[gtCavId].totalOverlap = 0.0;
			}
			numSetGridPoints = 0;
			for (wordId = 0; wordId < numBitWords; wordId++) { // overlap both grids and count all common grid points between the current MR cavity and the current GT cavity
				gtBitWordsGrid[wordId] &= mrBitWordsGrid[wordId];
				if (gtBitWordsGrid[wordId] == (uint64_t)0) continue;
				numSetGridPoints += CountBitsSet(gtBitWordsGrid[wordId]);
			}
			//sumOverlaps = (((double)numSetGridPoints)*gridPointVolume) / (gtCavities[gtCavId].totalVolume);
			sumOverlaps = (double)numSetGridPoints;
			gtCavities[gtCavId].totalOverlap += sumOverlaps;
			mrCavities[mrCavId].totalOverlap += sumOverlaps;
			//fprintf(outputFile, " %lf", sumOverlaps);
			fprintf(outputFile, " %lf", (sumOverlaps / gtCavities[gtCavId].totalVolume));
			memset(gtBitWordsGrid, 0, numBitWords * sizeof(uint64_t)); // reset grid of GT cavity to all zeroes
		} // loop for all GT cavities
		fprintf(outputFile, "\n");
		memset(mrBitWordsGrid, 0, numBitWords * sizeof(uint64_t)); // reset MR grid to all zeroes
	} // loop for all MR cavities
	fclose(outputFile);
	numTP = 0; // true positives (detected MR cavities overlapping true GT cavities)
	numFN = 0; // false negatives (true GT cavities not detected by any MR cavities)
	pctTP = 0.0;
	pctFN = 0.0;
	for (gtCavId = 0; gtCavId != numGTCavities; gtCavId++) {
		if (gtCavities[gtCavId].totalOverlap != 0.0) numTP++;
		else numFN++;
		pctTP += gtCavities[gtCavId].totalOverlap;
		pctFN += (gtCavities[gtCavId].totalVolume - gtCavities[gtCavId].totalOverlap);
	}
	pctTP = (pctTP / (double)totalGTCavitiesVolume);
	pctFN = (pctFN / (double)totalGTCavitiesVolume);
	numFP = 0; // false positives (detected MR cavities with no overlap with true GT cavities)
	pctTN = 0.0; // true negatives (in percentage) (if all the MR cavities coincided exactly with the GT cavities, the TNs outside the intersection would be 100%)
	pctFP = 0.0;
	for (mrCavId = 0; mrCavId < numMRCavities; mrCavId++) {
		if (mrCavities[mrCavId].totalOverlap == 0.0) numFP++;
		pctFP += (mrCavities[mrCavId].totalVolume - mrCavities[mrCavId].totalOverlap);
		pctTN += mrCavities[mrCavId].totalOverlap; // total overlap of all MR cavities
	}
	pctTN = (pctTN / (double)totalMRCavitiesVolume);
	pctFP = (pctFP / (double)totalMRCavitiesVolume);
	free(gtBitWordsGrid);
	free(mrBitWordsGrid);
	//fprintf(statsFile, "%s,%d,%d,%d,%d,%d,%.3lf\n", fileId, numGTCavities, numMRCavities, numTP, numFP, numFN, pctTN); // "PDB,#GT,#MR,TP,FP,FN,TN\n"
	fprintf(statsFile, "%s,%d,%d,%.3lf,%.3lf,%.3lf,%.3lf\n", fileId, numGTCavities, numMRCavities, pctTP, pctFP, pctFN, pctTN); // "PDB,#GT,#MR,TP,FP,FN,TN\n"
	free(fileId);
	if (createVis == 0) return 0;
	visualOutputFilename = CreateFilename(filename, outputDir, "-visualization.m", 5, 2);
	if ((visualOutputFile = fopen(visualOutputFilename, "w")) == NULL) {
		printf("\n> ERROR: Cannot create visual output file \"%s\"\n", visualOutputFilename);
		return (-1);
	}
	free(visualOutputFilename);
	fprintf(visualOutputFile, "(* ::Package:: *)\n\n");
	visualGridSize = 1.0f; // create visualization with a larger grid size
	createVis = 2;
	visualGridSize = (1.0f / 2.0f);
	minX -= (visualGridSize - gridSize);
	minY -= (visualGridSize - gridSize);
	minZ -= (visualGridSize - gridSize);
	gridSizeX = (unsigned int)ceilf((rangeX + 2 * (visualGridSize - gridSize)) / visualGridSize); // update new grid size in ranges
	gridSizeY = (unsigned int)ceilf((rangeY + 2 * (visualGridSize - gridSize)) / visualGridSize);
	gridSizeZ = (unsigned int)ceilf((rangeZ + 2 * (visualGridSize - gridSize)) / visualGridSize);
	minCoords[0] = minX;
	minCoords[1] = minY;
	minCoords[2] = minZ;
	gridSizes[0] = gridSizeX;
	gridSizes[1] = gridSizeY;
	gridSizes[2] = gridSizeZ;
	numBitWords = ((gridSizeX * gridSizeY * gridSizeZ) + 63) / 64; // resize grid arrays
	gtBitWordsGrid = (uint64_t *)calloc(numBitWords, sizeof(uint64_t));
	mrBitWordsGrid = (uint64_t *)calloc(numBitWords, sizeof(uint64_t));
	fprintf(visualOutputFile, "(* ::Input:: *)\n\n");
	for (mrCavId = 0; mrCavId < numMRCavities; mrCavId++) { // get a grid with all points from all MR cavities
		if (mrCavities[mrCavId].totalOverlap == 0.0) continue;
		LoadCavityIntoGrid(&(mrCavities[mrCavId]), mrBitWordsGrid, visualGridSize, gridSizes, minCoords); // the bits are "OR"-ed into the existing grid
	}
	for (gtCavId = 0; gtCavId != numGTCavities; gtCavId++) { // print the grid points of all GT cavities intersected by MR cavities
		LoadCavityIntoGrid(&(gtCavities[gtCavId]), gtBitWordsGrid, visualGridSize, gridSizes, minCoords);
		fprintf(visualOutputFile, "OverlappedGT[%d]={", (gtCavId + 1));
		numSetGridPoints = 0; // just to count the number of points in the intersection
		for (wordId = 0; wordId < numBitWords; wordId++) {
			gtBitWordsGrid[wordId] &= mrBitWordsGrid[wordId]; // get the common/shared points between the GT cavity and the MR cavities
			if (gtBitWordsGrid[wordId] == (uint64_t)0) continue;
			gridPoint = (wordId * 64);
			tempBitWord = gtBitWordsGrid[wordId];
			while (tempBitWord != 0ULL) { // get the coordinates of each one of the set grid points
				if (tempBitWord & 1ULL) { // check if rightmost bit is set
					divres = ldiv(gridPoint, gridSizeX); // gridPoint = (z*gridSizeY + y)*gridSizeX + x;
					x = (minX + visualGridSize * divres.rem); // convert from grid coordinates to real coordinates
					divres = ldiv(divres.quot, gridSizeY);
					y = (minY + visualGridSize * divres.rem);
					z = (minZ + visualGridSize * divres.quot);
					if (createVis == 1) fprintf(visualOutputFile, "%sPoint[{%.3f,%.3f,%.3f}]", ((numSetGridPoints == 0) ? "" : ","), x, y, z);
					else fprintf(visualOutputFile, "%sCuboid[{%.3f,%.3f,%.3f},{%.3f,%.3f,%.3f}]", ((numSetGridPoints == 0) ? "" : ","), x, y, z, (x + visualGridSize), (y + visualGridSize), (z + visualGridSize));
					numSetGridPoints++;
				}
				tempBitWord >>= 1; // next bit/point
				gridPoint++;
			}
		}
		fprintf(visualOutputFile, "};\n");
		memset(gtBitWordsGrid, 0, numBitWords * sizeof(uint64_t)); // now reset the grid of the GT cavity and get the non-overlapped points
		LoadCavityIntoGrid(&(gtCavities[gtCavId]), gtBitWordsGrid, visualGridSize, gridSizes, minCoords); // load the GT cavity again
		fprintf(visualOutputFile, "NotOverlappedGT[%d]={", (gtCavId + 1));
		numSetGridPoints = 0;
		for (wordId = 0; wordId < numBitWords; wordId++) {
			gtBitWordsGrid[wordId] &= ~(mrBitWordsGrid[wordId]); // get the GT points that were not overlapped
			if (gtBitWordsGrid[wordId] == (uint64_t)0) continue;
			gridPoint = (wordId * 64);
			tempBitWord = gtBitWordsGrid[wordId];
			while (tempBitWord != 0ULL) { // get the coordinates of each one of the set grid points
				if (tempBitWord & 1ULL) { // check if rightmost bit is set
					divres = ldiv(gridPoint, gridSizeX); // gridPoint = (z*gridSizeY + y)*gridSizeX + x;
					x = (minX + visualGridSize * divres.rem); // convert from grid coordinates to real coordinates
					divres = ldiv(divres.quot, gridSizeY);
					y = (minY + visualGridSize * divres.rem);
					z = (minZ + visualGridSize * divres.quot);
					if (createVis == 1) fprintf(visualOutputFile, "%sPoint[{%.3f,%.3f,%.3f}]", ((numSetGridPoints == 0) ? "" : ","), x, y, z);
					else fprintf(visualOutputFile, "%sCuboid[{%.3f,%.3f,%.3f},{%.3f,%.3f,%.3f}]", ((numSetGridPoints == 0) ? "" : ","), x, y, z, (x+visualGridSize), (y+visualGridSize), (z+visualGridSize));
					numSetGridPoints++;
				}
				tempBitWord >>= 1; // next bit/point
				gridPoint++;
			}
		}
		fprintf(visualOutputFile, "};\n");
		memset(gtBitWordsGrid, 0, numBitWords * sizeof(uint64_t)); // reset the grid for the next GT cavity
	}
	free(gtBitWordsGrid);
	free(mrBitWordsGrid);
	fprintf(visualOutputFile, "(* ::Input:: *)\n\n");
	fprintf(visualOutputFile, "Graphics3D[{\n");
	fprintf(visualOutputFile, "EdgeForm[],PointSize[%s],\n", (visualGridSize >= 1.0) ? "Large" : "Medium");
	//fprintf(visualOutputFile, "Black,Opacity[0.05],\n");
	//fprintf(visualOutputFile, "NotOverlappedMR,\n");
	for (gtCavId = 0; gtCavId != numGTCavities; gtCavId++) {
		fprintf(visualOutputFile, "Hue[%.2f],Opacity[%s],NotOverlappedGT[%d]", ((float)gtCavId / (float)numGTCavities), (createVis == 1)?("0.05"):("0.01"), (gtCavId + 1));
		if (gtCavities[gtCavId].totalOverlap != 0.0) fprintf(visualOutputFile, ",Opacity[1.0],OverlappedGT[%d]", (gtCavId + 1));
		if (gtCavId != (numGTCavities - 1)) fprintf(visualOutputFile, ",");
		fprintf(visualOutputFile, "\n");
	}
	fprintf(visualOutputFile, "},Boxed->True,Axes->True,AxesLabel->{\"X\",\"Y\",\"Z\"},LabelStyle->Directive[Bold],ImageSize->UpTo[2048]\n");
	fprintf(visualOutputFile, ",PlotLabel->\"#Cavities = %d : %d\nTP = %d ; FN = %d ; FP = %d\nTN = ~%.3lf\"\n", numGTCavities, numMRCavities, numTP, numFN, numFP, pctTN);
	fprintf(visualOutputFile, ",BaseStyle->{Antialiasing->True,RenderingOptions->{\"3DRenderingMethod\"->\"BSPTree\"}}\n"); // to better render transparent objects
	fprintf(visualOutputFile, "]\n");
	fprintf(visualOutputFile, "(* ::Input:: *)\n\n");
	fprintf(visualOutputFile, "Export[NotebookFileName[]<>\"_grid_%.2f_%s.png\",%%]\n", visualGridSize, (createVis == 1) ? "points" : "cubes");
	fclose(visualOutputFile);
	return 0;
}


// Code bellow for the old approach based on the mathematical formula
// for the overlap volume between two spheres
/*
// calculate 3D distance between the center of two spheres
float SpheresDistance(SphereCoords *s1, SphereCoords *s2) {
	float coordDist;
	float dist = 0.0;
	coordDist = ((s1->x) - (s2->x));
	dist += coordDist*coordDist;
	coordDist = ((s1->y) - (s2->y));
	dist += coordDist*coordDist;
	coordDist = ((s1->z) - (s2->z));
	dist += coordDist*coordDist;
	dist = sqrtf(dist);
	return dist;
}

int ProcessCavities(char *filename) {
	char *outputFilename, *visualOutputFilename;
	FILE *outputFile, *visualOutputFile;
	int numGTCavitySpheres, numMRCavitySpheres;
	int gtCavId, gtSphereId, mrCavId, mrSphereId;
	float distance;
	int numTP, numFP, numFN;
	double numTN, sumOverlaps;
	SphereCoords *gtSphere, *mrSphere;
	outputFilename = AddToFilename(filename, ".\\", "-output.txt", 1, 0);
	if ((outputFile = fopen(outputFilename, "w")) == NULL) {
		printf("\n> ERROR: Cannot create output file \"%s\"\n", outputFilename);
		return (-1);
	}
	free(outputFilename);
	visualOutputFilename = AddToFilename(filename, ".\\", "-visualization.m", 1, 0);
	if ((visualOutputFile = fopen(visualOutputFilename, "w")) == NULL) {
		printf("\n> ERROR: Cannot create visual output file \"%s\"\n", visualOutputFilename);
		return (-1);
	}
	free(visualOutputFilename);
	fprintf(visualOutputFile, "(* ::Package:: *)\n\n");
	fprintf(visualOutputFile, "(* ::Input:: *)\n\n");
	fprintf(outputFile, "%s", filename);
	for (gtCavId = 0; gtCavId < numGTCavities; gtCavId++) { // loop through all GT cavities
		numGTCavitySpheres = gtCavities[gtCavId].numSpheres;
		gtCavities[gtCavId].totalVolume = ((4.0 * M_PI * powf(1.0, 3.0)) / 3.0) * numGTCavitySpheres; // (4/3 * PI * r^3) * #spheres
		gtCavities[gtCavId].totalOverlap = 0.0;
		fprintf(visualOutputFile, "CavitiesGT[%d]={", (gtCavId + 1));
		for (gtSphereId = 0; gtSphereId < numGTCavitySpheres; gtSphereId++) {
			gtSphere = &(gtCavities[gtCavId].spheres[gtSphereId]);
			if (gtSphereId != 0) fprintf(visualOutputFile, ",");
			fprintf(visualOutputFile, "Ball[{%.3f,%.3f,%.3f},1.0]", (gtSphere->x), (gtSphere->y), (gtSphere->z));
		}
		fprintf(visualOutputFile, "};\n");
		fprintf(outputFile, " %d", gtCavId);
	}
	fprintf(outputFile, "\n");
	fprintf(visualOutputFile, "(* ::Input:: *)\n\n");
	for (mrCavId = 0; mrCavId < numMRCavities; mrCavId++) { // loop through all MR cavities
		numMRCavitySpheres = mrCavities[mrCavId].numSpheres;
		mrCavities[mrCavId].totalVolume = ((4.0 * M_PI * powf(1.0, 3.0)) / 3.0) * numMRCavitySpheres; // (4/3 * PI * r^3) * #spheres
		mrCavities[mrCavId].totalOverlap = 0.0;
		fprintf(visualOutputFile, "CavitiesMR[%d]={", (mrCavId + 1));
		for (mrSphereId = 0; mrSphereId < numMRCavitySpheres; mrSphereId++) {
			mrSphere = &(mrCavities[mrCavId].spheres[mrSphereId]);
			if (mrSphereId != 0) fprintf(visualOutputFile, ",");
			fprintf(visualOutputFile, "Ball[{%.3f,%.3f,%.3f},1.0]", (mrSphere->x), (mrSphere->y), (mrSphere->z));
		}
		fprintf(visualOutputFile, "};\n");
		fprintf(outputFile, "%d", mrCavId);
		for (gtCavId = 0; gtCavId < numGTCavities; gtCavId++) { // check overlap of this MR cavity with all GT cavities
			sumOverlaps = 0.0; // overlap between the current MR cavity and the current GT cavity
			numGTCavitySpheres = gtCavities[gtCavId].numSpheres;
			for (gtSphereId = 0; gtSphereId < numGTCavitySpheres; gtSphereId++) { // for each GT sphere, check overlap with all spheres of the current MR cavity
				gtSphere = &(gtCavities[gtCavId].spheres[gtSphereId]);
				for (mrSphereId = 0; mrSphereId < numMRCavitySpheres; mrSphereId++) {
					mrSphere = &(mrCavities[mrCavId].spheres[mrSphereId]);
					distance = SpheresDistance(mrSphere, gtSphere);
					if (distance < 2.0) { // intersecting spheres
						sumOverlaps += (M_PI * (4.0*1.0 + distance)*(2.0*1.0 - distance)*(2.0*1.0 - distance)) / 12.0; // 1/12 * PI * (4R+d) * (2R-d)^2
						(gtSphere->intersected) = 1.0; // mark GT sphere as intersected by MR sphere and vice-versa
						(mrSphere->intersected) = 1.0;
					}
				} // loop for MR spheres
			} // loop for GT spheres
			sumOverlaps = (sumOverlaps / gtCavities[gtCavId].totalVolume); // normalize over total GT cavity volume
			fprintf(outputFile, " %lf", sumOverlaps);
			gtCavities[gtCavId].totalOverlap += sumOverlaps; // update total overlap of this GT cavity (across all MR cavities), and vice-versa
			mrCavities[mrCavId].totalOverlap += sumOverlaps;
		} // loop for GT cavities
		fprintf(outputFile, "\n");
	} // loop for all MR cavities
	fclose(outputFile);
	numFP = 0; // false positives (detected MR cavities with no overlap with true GT cavities)
	for (mrCavId = 0; mrCavId < numMRCavities; mrCavId++) {
		if (mrCavities[mrCavId].totalOverlap == 0.0) numFP++;
	}
	numTP = 0; // true positives (detected MR cavities overlapping true GT cavities)
	numFN = 0; // false negatives (true GT cavities not detected by any MR cavities)
	sumOverlaps = 0.0; // total overlap of all GT cavities
	for (gtCavId = 0; gtCavId != numGTCavities; gtCavId++) {
		if (gtCavities[gtCavId].totalOverlap != 0.0) numTP++;
		else numFN++;
		sumOverlaps += gtCavities[gtCavId].totalOverlap;
	}
	numTN = (sumOverlaps / (double)numGTCavities); // true negatives (if all the MR cavities coincided exactly with the GT cavities, the TNs outside the intersection would be 100%)
	fprintf(visualOutputFile, "(* ::Input:: *)\n\n");
	for (gtCavId = 0; gtCavId != numGTCavities; gtCavId++) { // print all GT spheres intersected by MR spheres
		fprintf(visualOutputFile, "Overlapped[%d]={", (gtCavId + 1));
		numGTCavitySpheres = gtCavities[gtCavId].numSpheres;
		sumOverlaps = 0.0; // just to count the number of overlapped spheres
		for (gtSphereId = 0; gtSphereId < numGTCavitySpheres; gtSphereId++) {
			gtSphere = &(gtCavities[gtCavId].spheres[gtSphereId]);
			if ((gtSphere->intersected) == 0.0) continue;
			fprintf(visualOutputFile, "%sBall[{%.3f,%.3f,%.3f},1.0]", ((sumOverlaps == 0.0) ? "" : ","), (gtSphere->x), (gtSphere->y), (gtSphere->z));
			sumOverlaps++;
		}
		fprintf(visualOutputFile, "};\n");
	}
	fprintf(visualOutputFile, "(* ::Input:: *)\n\n");
	fprintf(visualOutputFile, "Style[\n");
	fprintf(visualOutputFile, "Graphics3D[{\n");
	fprintf(visualOutputFile, "Black,Opacity[0.05],\n");
	for (mrCavId = 0; mrCavId != numMRCavities; mrCavId++) {
		fprintf(visualOutputFile, "CavitiesMR[%d],\n", (mrCavId + 1));
	}
	for (gtCavId = 0; gtCavId != numGTCavities; gtCavId++) {
		fprintf(visualOutputFile, "Hue[%.2f],Opacity[0.05],CavitiesGT[%d]", ((float)gtCavId / (float)numGTCavities), (gtCavId + 1));
		if (gtCavities[gtCavId].totalOverlap != 0.0) fprintf(visualOutputFile, ",Opacity[1.0],Overlapped[%d]", (gtCavId + 1));
		if (gtCavId != (numGTCavities - 1)) fprintf(visualOutputFile, ",");
		fprintf(visualOutputFile, "\n");
	}
	fprintf(visualOutputFile, "},Boxed->True,Axes->True,AxesLabel->{x,y,z},LabelStyle->Directive[Bold],ImageSize->UpTo[1024]");
	fprintf(visualOutputFile, ",PlotLabel->\"#Cavities = %d : %d\nTP = %d ; FN = %d ; FP = %d\nTN = ~%.3lf\"", numGTCavities, numMRCavities, numTP, numFN, numFP, numTN);
	fprintf(visualOutputFile, "]\n");
	//fprintf(visualOutputFile, ",RenderingOptions->{\"3DRenderingMethod\"->\"BSPTree\"}");
	fprintf(visualOutputFile, "]\n");
	fclose(visualOutputFile);
	return 0;
}
*/

int main(int argc, char *argv[]) {
	int numFiles, numErrors;
	char *gtDirName, *mrDirName;
	char *gtFilename, *mrFilename;
	char *statsFilename, *errorsFilename;
	char *filename;
	int withVis;
#ifdef _WIN32
	WIN32_FIND_DATA FindFileData;
	HANDLE hFind;
	int i;
#else
	DIR *dirStream;
	struct dirent *dirEntry;
	char *tempExtension;
#endif

	printf("[ CavitiesOverlap v0.5.2 ]\n");
	if (argc<3) {
		printf("\nUsage:\n");
		printf("\t%s <Ground_Truth_Directory> <Method_Results_Directory> (<Output_Directory>)\n", argv[0]);
		printf("\nHelp:\n");
		printf("\t<Ground_Truth_Directory>:\tdirectory containing .CSV files (required)\n");
		printf("\t<Method_Results_Directory>:\tdirectory containing .TXT files (required)\n");
		printf("\t<Output_Directory>:\t\tdirectory to save output (default is MR directory) (optional)\n");
		printf("\nOptions:\n");
		printf("\t-v\tGenerate visualizations (in 'Mathematica Notebook' format)\n");
		printf("\n");
		return (-1);
	}

#ifdef _WIN32
	mrDirName = CreateFilename(argv[2], NULL, "*", 3, 0); // add wildcard to the MR directory
	if ((hFind = FindFirstFile(mrDirName, &FindFileData)) == INVALID_HANDLE_VALUE) {
		printf("> ERROR: Failed to locate method results file(s) in \"%s\"\n", mrDirName);
		return (-1);
	}
	FindClose(hFind);
	gtDirName = CreateFilename(argv[1], NULL, "*.csv", 3, 0); // add wildcard to the GT directory
	if ((hFind = FindFirstFile(gtDirName, &FindFileData)) == INVALID_HANDLE_VALUE) { // set GT dir after MR dir so 'hFind' var enters the main loop with this dir loaded
		printf("> ERROR: Failed to locate ground truth file(s) in \"%s\"\n", gtDirName);
		return (-1);
	}
	for (i = 0; mrDirName[i] != '*'; i++); // remove the '*' from the directories' names
	mrDirName[i] = '\0';
	for (i = 0; gtDirName[i] != '*'; i++);
	gtDirName[i] = '\0';
#else
	mrDirName = CreateFilename(argv[2], NULL, NULL, 3, 0);
	if ((dirStream = opendir(mrDirName)) == NULL) {
		printf("> ERROR: Failed to locate method results directory \"%s\"\n", mrDirName);
		return (-1);
	}
	closedir(dirStream);
	gtDirName = CreateFilename(argv[1], NULL, NULL, 3, 0);
	if ((dirStream = opendir(gtDirName)) == NULL) {
		printf("> ERROR: Failed to locate ground truth directory \"%s\"\n", gtDirName);
		return (-1);
	}
	tempExtension = NULL;
#endif
	
	/*
	withVis = 0;
	outputDir = NULL;
	for (i = 3; i < argc; i++) { // process command line arguments
		if (argv[i][0] != '-') continue;
		if (argv[i][1] == 'v' || argv[i][1] == 'V') withVis = 1;
		else if (argv[i][1] == 'o' || argv[i][1] == 'O') {
			if ((++i) == argc) break;
			outputDir = CreateFilename(argv[i], NULL, NULL, 3, 0); // add a slash to the output path if needed
		}
	}
	*/
	if (argc >= 5 && argv[4][0] == '-' && (argv[4][1] == 'v' || argv[4][1] == 'V')) withVis = 1; // process 4th command line argument
	else withVis = 0;
	if (argc >= 4 && argv[3][0] != '-') outputDir = CreateFilename(argv[3], NULL, NULL, 3, 0); // add a slash to the output path if needed
	else outputDir = CreateFilename(mrDirName, NULL, NULL, 3, 0); // if no output path was given, set it to the MR directory
	statsFilename = CreateFilename(outputDir, NULL, FILENAME_STATS, 0, 0);
	errorsFilename = CreateFilename(outputDir, NULL, FILENAME_ERRORS, 0, 0);

	if ((statsFile = fopen(statsFilename, "w")) == NULL) {
		printf("> ERROR: Cannot create statistics file \"%s\"\n", statsFilename);
		return (-1);
	}
	fprintf(statsFile, "PDB,#GT,#MR,TP,FP,FN,TN\n");

	numFiles = 0;
	numErrors = 0;
	while (1) { // process all pairs of files
#ifdef _WIN32
		if (FindFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) continue;
		filename = FindFileData.cFileName;
#else
		if ((dirEntry = readdir(dirStream)) == NULL) break;
		if ((dirEntry->d_type) == DT_DIR) continue;
		filename = (dirEntry->d_name);
		if (tempExtension != NULL) free(tempExtension);
		tempExtension = CreateFilename(filename, NULL, NULL, 4, 0);
		if (strcmp(tempExtension, ".csv") != 0) continue;
#endif
		//if (strcmp(filename, FILENAME_STATS) == 0) continue; // on Windows, it will loop endlessly without fetching another file
		gtFilename = CreateFilename(filename, gtDirName, NULL, 0, 0);
		mrFilename = CreateFilename(filename, mrDirName, ".txt", 1, 0);
		numFiles++;
		printf("> Processing files <%s> <%s> ... ", gtFilename, mrFilename);
		//printf("> Processing files ");
		fflush(stdout);
		if ((numGTCavities = LoadCavitiesFile(gtFilename, &(gtCavities))) > 0 && (numMRCavities = LoadCavitiesFile(mrFilename, &(mrCavities))) > 0) {
			//printf("... ");
			ProcessCavitiesWithGrid(GRID_SIZE, mrFilename, withVis); // each file pair is processed here
			printf("\n");
		}
		else { // error loading the cavities of one of the two files
			numErrors++;
			if (errorsFile == NULL) {
				if ((errorsFile = fopen(errorsFilename, "w")) == NULL) {
					printf("> ERROR: Cannot create errors file \"%s\"\n", errorsFilename);
					return (-1);
				}
				fprintf(errorsFile, "File,Line\n");
			}
			if (numGTCavities <= 0) fprintf(errorsFile, "%s,%d\n", gtFilename, (-numGTCavities));
			else fprintf(errorsFile, "%s,%d\n", mrFilename, (-numMRCavities));
		}
		FreeCavityArrays();
		free(gtFilename);
		free(mrFilename);
#ifdef _WIN32
		if (FindNextFile(hFind, &FindFileData) == 0) break;
#endif
	}
	printf("> %d total file pairs processed", numFiles);
	if (numErrors != 0) printf(" (%d errors)", numErrors);
	printf("\n");
	printf("> Saving statistics to \"%s\"\n", statsFilename);
	fclose(statsFile);
	if (errorsFile != NULL) {
		printf("> Saving error log to \"%s\"\n", errorsFilename);
		fclose(errorsFile);
	}

#ifdef _WIN32
	FindClose(hFind);
#else
	closedir(dirStream);
	if (tempExtension != NULL) free(tempExtension);
#endif

	free(mrDirName);
	free(gtDirName);
	free(outputDir);
	free(statsFilename);
	free(errorsFilename);
	
	printf("> Done!\n");
#if PAUSE_AT_EXIT == 1
	getchar();
#endif
	return 0;
}
