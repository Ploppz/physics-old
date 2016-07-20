#include "FontTexture.h"
#include "skyline-bl.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "edtaa3func.h"

#include <ft2build.h>
#ifdef FT_FREETYPE_H
#include FT_FREETYPE_H
// #else
// #warning "Not including."
#endif


#include <assert.h>
#include <string>
#include <iostream>
#include <ios>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstdio>


void insertImage(unsigned char *dest, int destW, int destH, unsigned char* src, int srcW, int srcH, int srcX, int srcY);
FT_Bitmap copyBitmap(FT_Bitmap &bm);
float *computeDF(unsigned char *image, int width, int height);
void renderImage(float *image, int width, int height);
unsigned char *invertImage(unsigned char *image, int width, int height);
FT_Bitmap renderBitmap(FT_Bitmap bitmap);

FontTexture::FontTexture()
	: rawAtlas(0), atlas{}, glyphs{}
{
	
}

Glyph &FontTexture::getGlyph(unsigned char ch)
{
	return glyphs[ch];
}

template <typename T>
void writeArrayToFile(std::ofstream &f, std::vector<T> &array)
{
	f.write((char *)array.data(), array.size() * sizeof(T));
}
template <typename T>
void readArrayFromFile(std::ifstream &f, std::vector<T> &array)
{
	// Length is here number of bytes
	int startPos = f.tellg();
	f.seekg(0, std::ios::end);
	int length = int(f.tellg()) - startPos;
	f.seekg(startPos, std::ios::beg);
	//
	array.resize(length / sizeof(T));
	f.read((char *)array.data(), length);
}

void FontTexture::useExistingAtlas(std::string atlaspath, std::string metadatapath)
{
	std::cout << "Call to useExistingAtlas" << std::endl;

	// Open files
	std::ifstream atlas_f;
	atlas_f.open(atlaspath, std::ios::binary | std::ios::in);
	if (!atlas_f)  {std::cerr << "Could not open file " << atlaspath << std::endl; exit(1);}
	std::ifstream meta_f;
	meta_f.open(metadatapath, std::ios::binary | std::ios::in);

	if (!meta_f)   {std::cerr << "Could not open file " << metadatapath << std::endl; exit(1);}

	// First read width and height
	atlas_f.read((char *)(&atlasWidth), sizeof(atlasWidth));
	atlas_f.read((char *)(&atlasHeight), sizeof(atlasHeight));
	readArrayFromFile(atlas_f, atlas);
	readArrayFromFile(meta_f, glyphs);

	atlas_f.close();
	meta_f.close();

}
void FontTexture::writeAtlasToFile(std::string outatlas, std::string outmeta)
{
	std::cout << "Call to writeAtlasToFile" << std::endl;

	// Open files
	std::ofstream atlas_f;
	atlas_f.open(outatlas, std::ios::binary | std::ios::out);
	if (!atlas_f)  {std::cerr << "Could not open file " << outatlas << std::endl; exit(1);}
	std::ofstream meta_f;
	meta_f.open(outmeta, std::ios::binary | std::ios::out);
	if (!meta_f)   {std::cerr << "Could not open file " << outmeta << std::endl; exit(1);}

	// First write width and height to file
	atlas_f.write((char *)(&atlasWidth), sizeof(atlasWidth));
	atlas_f.write((char *)(&atlasHeight), sizeof(atlasHeight));
	writeArrayToFile(atlas_f, atlas);
	writeArrayToFile(meta_f, glyphs);

	atlas_f.close();
	meta_f.close();
}


// Only used for algorithm - as a key. Reason: To access the bitmap of the glyph, and update the texture coordinates.
struct GlyphHolder
{
	Glyph *glyph;
	FT_Bitmap bitmap;
};


void FontTexture::generateAtlas(std::string fontpath)
{
	generateAtlas(fontpath, 1, NORMAL);
}


/*
	Generate image and metadata, store in RAM. Make a function that saves these two things to storage.
*/
void FontTexture::generateAtlas(std::string fontpath, int resolution, AtlasType type)
{
	FT_Library FTlib;
	FT_Face face;
	int error = FT_Init_FreeType( &FTlib);
	assert( !error );

	this->resolution = resolution;
	this->type = type;

	// Load font
	error = FT_New_Face(FTlib, fontpath.c_str(), 0, &face);
	if (error == FT_Err_Unknown_File_Format) {
		std::cout << fontpath << ": Unsupported font format." << std::endl;
	} else if (error) {
		std::cout << "Error opening file " << fontpath << ". Broken?" << std::endl;
		exit(1);
	}
	if (resolution != 1) {
		if(FT_Set_Pixel_Sizes(face, 0, resolution)) assert(!"Wrong size");
	}


// 1. Render glyphs individually and save into vectors.

	glyphs = std::vector<Glyph>(255);
	std::vector<SBL::Size<GlyphHolder>> glyphBoxes; // Used only to send into Skyline algorithm

	// Render (FOR NOW) glyphs individually, from 32 to 254
	SBL::Size<GlyphHolder> newBox;
	Glyph newGlyph;
	FT_UInt charIndex, charIndex2;
	FT_Vector kerning;
	int i = 0;
	for (unsigned char c = 32; c <= 254; c ++, i ++)
	{
		// If newGlyph.charcode is still 0 in the end, then this entry is trash
		newGlyph = Glyph {};

		charIndex = FT_Get_Char_Index( face, c );
		if (!charIndex) continue;

		if (!FT_Load_Char(face, c, 0)) {
			if (!FT_Render_Glyph(face->glyph, FT_RENDER_MODE_NORMAL)) {
				// Data & metrics are now in face->glyph (type FT_GlyphSlot)
				newGlyph.twidth = 	face->glyph->metrics.width / 64.0f;
				newGlyph.theight = 	face->glyph->metrics.height / 64.0f;
				newGlyph.width = static_cast<float>(newGlyph.twidth) / resolution;
				newGlyph.height = static_cast<float>(newGlyph.theight) / resolution;
				newGlyph.xadvance = face->glyph->metrics.horiAdvance / 64.0f / resolution;
				// newGlyph.yadvance = face->glyph->advance.y / 64.0f / resolution;
				// newGlyph.xoffset = face->glyph->metrics.horiBearingX / 64.0f / resolution;
				// newGlyph.yoffset = face->glyph->metrics.horiBearingY / 64.0f / resolution;
				newGlyph.xoffset = face->glyph->bitmap_left / static_cast<float>(resolution);
				newGlyph.yoffset = -face->glyph->bitmap_top/ static_cast<float>(resolution);
				newGlyph.charcode = c;
				
				
				newBox.width = newGlyph.twidth + 10;
				newBox.height = newGlyph.theight + 10;

				newBox.key.glyph = &glyphs[c];
				newBox.key.bitmap = renderBitmap(face->glyph->bitmap); // Need deep-copy because it gets destructed.

				glyphBoxes.push_back(newBox);
			} else {
				std::cerr << "Could not render glyph " << c << " (" << static_cast<int>(c) << ")." << std::endl;
			}
		} else {
			std::cerr << "Could not load char " << c << " (" << static_cast<int>(c) << ")." << std::endl;
		}
		glyphs[c] = newGlyph;
	}

	// Kerning
	FT_UInt leftIndex, rightIndex;
	for (unsigned char left = 0; left <= 254; left ++)
	{
		leftIndex = FT_Get_Char_Index(face, left);
		for (unsigned char right = 0; right <= 254; right ++)
		{
			rightIndex = FT_Get_Char_Index(face, right);
			FT_Get_Kerning(face, leftIndex, rightIndex, FT_KERNING_DEFAULT, &kerning);
			// NOTE: Kerning of a glyph is with respect to preceeding glyphs when rendering
			glyphs[right].kerning[left] = kerning.x / 64.0f / resolution;
		}
	}

// 2. Use some of the previous results to put all the glyph boxes into one big box.

	// Send the font "boxes" (width, height, glyph ptr) to Skyline algorithm to pack into rectangle.
	double atlasWidthd, atlasHeightd;
	SBL::SkylineBL_alg<GlyphHolder> skylineAlg(glyphBoxes);
	// SBL::Area describes where each box resides (x, y, width, height, glyph ptr
	std::vector<SBL::Area<GlyphHolder>> virtualAtlas = skylineAlg.compute(atlasWidthd, atlasHeightd);
	atlasWidth = static_cast<int>(atlasWidthd);
	atlasHeight = static_cast<int>(atlasHeightd);

	// Note: Coordinates are with origin bottom-left


// 3. Render the atlas to an image
	int length = atlasWidth * atlasHeight;
	rawAtlas = new unsigned char[length];

	for (auto it = virtualAtlas.begin(); it != virtualAtlas.end(); it ++)
	{
		
		if (!(it->key.bitmap.width == 0 || it->key.bitmap.rows == 0)) {
			insertImage(rawAtlas, atlasWidth, atlasHeight, // dest
					it->key.bitmap.buffer, it->key.bitmap.width, it->key.bitmap.rows, // src
					it->x, it->y); // position

		}
		// Update glyphs with texture coordinates (from the virtual atlas)
		it->key.glyph->u = it->x;
		it->key.glyph->v = it->y;
	}

	stbi_write_png("atlas.png", atlasWidth, atlasHeight, 1, rawAtlas, atlasWidth);
	

// 4. Apply distance transform
	
	if (type == SDF) {
		float *atlas_data = computeDF(rawAtlas, atlasWidth, atlasHeight);
		atlas = std::vector<float>(atlas_data, atlas_data + length); // Eventually use assign()
	} else {
		atlas = std::vector<float>(rawAtlas, rawAtlas + length);
	}
}

// Copies bitmap, and if each pixel is one bit, converts to one char per pixel
FT_Bitmap renderBitmap(FT_Bitmap bitmap)
{
	int size = bitmap.width * bitmap.rows;
	FT_Bitmap result = bitmap;
	result.buffer = new unsigned char[size];
	if (bitmap.pixel_mode == ft_pixel_mode_grays) {
		// Passthrough
		std::copy(bitmap.buffer, bitmap.buffer + size, result.buffer);
	} else if (bitmap.pixel_mode == ft_pixel_mode_mono) {
		// Convert to one char per pixel
		int byte = 0, bit = 7;
		int bitData;
		int treshold = 0;
		for (int i = 0; i < size; i ++) {
			if (bit < treshold) {
				bit = 7;
				byte ++;
			}
			// Change treshold if byte is right-most in the picture
			if (byte % bitmap.pitch == bitmap.pitch - 1) {
				treshold = bitmap.pitch * 8 - bitmap.width;
			} else treshold = 0;
			
			bitData = bitmap.buffer[byte] & (1 << bit);
			result.buffer[i] = bitData * 255;

			bit --;
		}
	}
	return result;
}

void insertImage(unsigned char *dest, int destW, int destH, unsigned char* src, int srcW, int srcH, int destX, int destY)
{
	assert(destX >= 0 && destY >= 0 && (destX + srcW) <= destW && (destY + srcH) <= destH);

	// For each row, copy.
	unsigned char *destination;
	unsigned char *sourceBase;
	int srcRow, destRow;
	for (destRow = destY, srcRow = 0; destRow < destH && srcRow < srcH; destRow ++, srcRow ++)
	{
		destination = dest + destX + destRow * destW;
		sourceBase = src + srcRow * srcW;
		std::copy(sourceBase, sourceBase + srcW, destination);	
	}
}


// Deep copy
FT_Bitmap copyBitmap(FT_Bitmap &bm)
{
	
	int size = bm.width * bm.rows;
	FT_Bitmap result = bm;
	result.buffer = new unsigned char[size];
	std::copy(bm.buffer, bm.buffer + size, result.buffer);

	return result;
}

FontTexture::~FontTexture()
{
	if (rawAtlas != 0) delete[] rawAtlas;
}


double *applyEdtaa3(unsigned char *image, int width, int height)
{
	int length = width * height;

	double *doubleData = new double[length];
	// Convert to double [0, 1]
	for (int i = 0; i < length; i ++)
	{
		doubleData[i] = static_cast<double>(image[i]) / 255.0f;
	}
	double *gx = new double[length];
	double *gy = new double[length];
	edtaa3::computegradient(doubleData, width, height, gx, gy);
	short *distx = new short[length];
	short *disty = new short[length];
	double *dist = new double[length];
	edtaa3::edtaa3(doubleData, gx, gy, width, height, distx, disty, dist);
	delete[] doubleData;
	delete[] distx;
	delete[] disty;
	delete[] gx;
	delete[] gy;
	// Make all pixels >= 0
	for (int i = 0; i < length; i ++)
	{
		if (dist[i] < 0) dist[i] = 0;
	}
	return dist;
}

float *toFloat(double *image, int size)
{
	float *ret = new float[size];
	for (int i = 0; i < size; i ++)
	{
		ret[i] = static_cast<float>(image[i]);
	}
	return ret;
}

// Compute distance field
float *computeDF(unsigned char *image, int width, int height)
{
	int length = width * height;
	double *dist1 = applyEdtaa3(image, width, height);
	unsigned char *invertedImage = invertImage(image, width, height);
	double *dist2 = applyEdtaa3(invertedImage, width, height);
	
	// "merge" the two distance fields: dist = dist1 - dist2
	double *dist = new double[length];
	for (int i = 0; i < length; i ++)
	{
		 dist[i] = dist1[i] - dist2[i];
	}

	float *floatData = toFloat(dist, width*height);

	delete[] dist1;
	delete[] dist2;
	delete[] dist;

	return floatData;
}
/* Allocates a new image rather than editing existing one. */
unsigned char *invertImage(unsigned char *image, int width, int height)
{
	int length = width * height;
	unsigned char *newImage = new unsigned char[length];
	for (int i = 0; i < length; i ++)
	{
		newImage[i] = 255 - image[i];
	}
	return newImage;
}
