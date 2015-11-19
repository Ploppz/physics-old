#pragma once

#include <ft2build.h>
#ifdef FT_FREETYPE_H
#include FT_FREETYPE_H
#else
#warning "Not including."
#endif

#include <string>
#include <vector>

struct Glyph
{
	Glyph() :charcode{}, width{}, height{}, xoffset{}, yoffset{}, u{}, v{}, twidth{}, theight{} {}
	wchar_t charcode;
	// Normalized:
	float xoffset, yoffset;
	float width, height;
	float xadvance, yadvance;
	//
	float kerning[255];
	
// Top left texture coordinate
	int u, v;
	// width and height in the texture
	int twidth, theight;
};

enum AtlasType {
	SDF, NORMAL
};

class FontTexture
{
private:

	int atlasWidth, atlasHeight;
	AtlasType type;
	int resolution = 16;

public:
	FontTexture();
	FontTexture(const FontTexture&) = delete;
	FontTexture(FontTexture &&) = delete;
	FontTexture& operator= (const FontTexture&) = delete;
	FontTexture& operator= (FontTexture&&) = delete;


	unsigned char *rawAtlas;
	std::vector<float> atlas;
	std::vector<Glyph> glyphs;

	Glyph &getGlyph(unsigned char ch);

	

	int getAtlasWidth() { return atlasWidth; }
	int getAtlasHeight() { return atlasHeight; }

	
	// If resolution == 1, doesn't rescale. Resolution is just the pixel size.
	void generateAtlas(std::string fontpath, int resolution, AtlasType type);
	// Use default font size, and dont' use SDF:
	void generateAtlas(std::string fontpath);
	void writeAtlasToFile(std::string outatlas, std::string outmeta);

	void useExistingAtlas(std::string atlaspath, std::string metadatapath);


	~FontTexture();
};


