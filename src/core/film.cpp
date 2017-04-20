
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// core/film.cpp*
#include "film.h"
#include "paramset.h"
#include "imageio.h"
#include "stats.h"

namespace pbrt {

STAT_MEMORY_COUNTER("Memory/Film pixels", filmPixelMemory);

// Film Method Definitions
Film::Film(const Point2i &resolution, const Bounds2f &cropWindow,
           std::unique_ptr<Filter> filt, Float diagonal,
           const std::string &filename, Float scale, Float maxSampleLuminance)
    : fullResolution(resolution),
      diagonal(diagonal * .001),
      filter(std::move(filt)),
      filename(filename),
      scale(scale),
      maxSampleLuminance(maxSampleLuminance) {
    // Compute film image bounds
    croppedPixelBounds =
        Bounds2i(Point2i(std::ceil(fullResolution.x * cropWindow.pMin.x),
                         std::ceil(fullResolution.y * cropWindow.pMin.y)),
                 Point2i(std::ceil(fullResolution.x * cropWindow.pMax.x),
                         std::ceil(fullResolution.y * cropWindow.pMax.y)));
    LOG(INFO) << "Created film with full resolution " << resolution <<
        ". Crop window of " << cropWindow << " -> croppedPixelBounds " <<
        croppedPixelBounds;

    // Allocate film image storage
    pixels = std::unique_ptr<Pixel[]>(new Pixel[croppedPixelBounds.Area()]);
    filmPixelMemory += croppedPixelBounds.Area() * sizeof(Pixel);

    // Precompute filter weight table
    int offset = 0;
    for (int y = 0; y < filterTableWidth; ++y) {
        for (int x = 0; x < filterTableWidth; ++x, ++offset) {
            Point2f p;
            p.x = (x + 0.5f) * filter->radius.x / filterTableWidth;
            p.y = (y + 0.5f) * filter->radius.y / filterTableWidth;
            filterTable[offset] = filter->Evaluate(p);
        }
    }
	if(PbrtOptions.gui)
	{
		int xPixelCount = fullResolution.x;
		int yPixelCount = fullResolution.y;
		SDL_Init(SDL_INIT_VIDEO);
		sdl_window = SDL_CreateWindow("My Game Window",
				SDL_WINDOWPOS_UNDEFINED,
				SDL_WINDOWPOS_UNDEFINED,
				xPixelCount, yPixelCount,
				SDL_WINDOW_OPENGL);
		if (sdl_window == NULL)
			printf("Unable to create window to display image");
		sdl_renderer = SDL_CreateRenderer(sdl_window, -1, 0);
		Uint32 rmask, gmask, bmask, amask;

		/* SDL interprets each pixel as a 32-bit number, so our masks must depend
		   on the endianness (byte order) of the machine */
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
		rmask = 0xff000000;
		gmask = 0x00ff0000;
		bmask = 0x0000ff00;
		amask = 0x000000ff;
#else
		rmask = 0x000000ff;
		gmask = 0x0000ff00;
		bmask = 0x00ff0000;
		amask = 0xff000000;
#endif

		/* Create a 32-bit surface with the bytes of each pixel in R,G,B,A order,
		   as expected by OpenGL for textures */
		//sdl_surface = SDL_CreateRGBSurface(0, xPixelCount, yPixelCount, 32, rmask, gmask, bmask, amask);
		//Using zeros for the RGB masks sets a default value, based on the depth. (e.g. SDL_CreateRGBSurface(0,w,h,32,0,0,0,0);) 
		//However, using zero for the Amask results in an Amask of 0
		sdl_surface = SDL_CreateRGBSurface(0, xPixelCount, yPixelCount, 32, 0, 0, 0, 0);
		if (sdl_surface == NULL) 
		{
			SDL_Log("SDL_CreateRGBSurface() failed: %s", SDL_GetError());
			exit(1);
		}
		sdl_texture = SDL_CreateTexture(sdl_renderer,
				SDL_PIXELFORMAT_ARGB8888,
				SDL_TEXTUREACCESS_STREAMING,
				xPixelCount, yPixelCount);
		if (SDL_LockSurface(sdl_surface) == -1)
		{
			printf("Unable to lock surface for image display");
		}
		else 
		{
			for (int y = 0; y < yPixelCount; ++y) 
			{
				for (int x = 0; x < xPixelCount; ++x) 
				{
					u_int *bufp = (u_int *)sdl_surface->pixels + y*sdl_surface->pitch/4 + x;
					*bufp = SDL_MapRGB(sdl_surface->format, 64, 64, 64);
				}
			}
			SDL_UnlockSurface(sdl_surface);
			SDL_SetRenderDrawColor(sdl_renderer, 255, 0, 0, 255);
			SDL_RenderClear(sdl_renderer);
			SDL_RenderPresent(sdl_renderer);// old: SDL_UpdateRect(sdlWindow, 0, 0, xPixelCount, yPixelCount);
			SDL_Delay(2000);
			//display frame
			SDL_SetRenderDrawColor(sdl_renderer, 0, 255, 0, 255);
			SDL_RenderClear(sdl_renderer);
			SDL_RenderPresent(sdl_renderer);// old: SDL_UpdateRect(sdlWindow, 0, 0, xPixelCount, yPixelCount);
			SDL_Delay(2000);
			//display frame
			SDL_SetRenderDrawColor(sdl_renderer, 0, 0, 255, 255);
			SDL_RenderClear(sdl_renderer);
			SDL_RenderPresent(sdl_renderer);// old: SDL_UpdateRect(sdlWindow, 0, 0, xPixelCount, yPixelCount);
			//display frame
			SDL_SetRenderDrawColor(sdl_renderer, 0, 0, 0, 255);
			SDL_RenderClear(sdl_renderer);
			SDL_RenderPresent(sdl_renderer);// old: SDL_UpdateRect(sdlWindow, 0, 0, xPixelCount, yPixelCount);
			SDL_Delay(2000);
			//display frame
		}
	}
	  }

Bounds2i Film::GetSampleBounds() const {
	Bounds2f floatBounds(Floor(Point2f(croppedPixelBounds.pMin) +
				Vector2f(0.5f, 0.5f) - filter->radius),
			Ceil(Point2f(croppedPixelBounds.pMax) -
				Vector2f(0.5f, 0.5f) + filter->radius));
	return (Bounds2i)floatBounds;
}

Bounds2f Film::GetPhysicalExtent() const {
	Float aspect = (Float)fullResolution.y / (Float)fullResolution.x;
    Float x = std::sqrt(diagonal * diagonal / (1 + aspect * aspect));
    Float y = aspect * x;
    return Bounds2f(Point2f(-x / 2, -y / 2), Point2f(x / 2, y / 2));
}

std::unique_ptr<FilmTile> Film::GetFilmTile(const Bounds2i &sampleBounds) {
    // Bound image pixels that samples in _sampleBounds_ contribute to
    Vector2f halfPixel = Vector2f(0.5f, 0.5f);
    Bounds2f floatBounds = (Bounds2f)sampleBounds;
    Point2i p0 = (Point2i)Ceil(floatBounds.pMin - halfPixel - filter->radius);
    Point2i p1 = (Point2i)Floor(floatBounds.pMax - halfPixel + filter->radius) +
                 Point2i(1, 1);
    Bounds2i tilePixelBounds = Intersect(Bounds2i(p0, p1), croppedPixelBounds);
    return std::unique_ptr<FilmTile>(new FilmTile(
        tilePixelBounds, filter->radius, filterTable, filterTableWidth,
        maxSampleLuminance));
}

void Film::Clear() {
    for (Point2i p : croppedPixelBounds) {
        Pixel &pixel = GetPixel(p);
        for (int c = 0; c < 3; ++c)
            pixel.splatXYZ[c] = pixel.xyz[c] = 0;
        pixel.filterWeightSum = 0;
    }
}

void Film::MergeFilmTile(std::unique_ptr<FilmTile> tile) {
    ProfilePhase p(Prof::MergeFilmTile);
    VLOG(1) << "Merging film tile " << tile->pixelBounds;
    std::lock_guard<std::mutex> lock(mutex);
    for (Point2i pixel : tile->GetPixelBounds()) {
        // Merge _pixel_ into _Film::pixels_
        const FilmTilePixel &tilePixel = tile->GetPixel(pixel);
        Pixel &mergePixel = GetPixel(pixel);
        Float xyz[3];
        tilePixel.contribSum.ToXYZ(xyz);
        for (int i = 0; i < 3; ++i) mergePixel.xyz[i] += xyz[i];
        mergePixel.filterWeightSum += tilePixel.filterWeightSum;
    }
}

void Film::SetImage(const Spectrum *img) const {
    int nPixels = croppedPixelBounds.Area();
    for (int i = 0; i < nPixels; ++i) {
        Pixel &p = pixels[i];
        img[i].ToXYZ(p.xyz);
        p.filterWeightSum = 1;
        p.splatXYZ[0] = p.splatXYZ[1] = p.splatXYZ[2] = 0;
    }
}

void Film::AddSplat(const Point2f &p, Spectrum v) {
    ProfilePhase pp(Prof::SplatFilm);

    if (v.HasNaNs()) {
        LOG(ERROR) << StringPrintf("Ignoring splatted spectrum with NaN values "
                                   "at (%f, %f)", p.x, p.y);
        return;
    } else if (v.y() < 0.) {
        LOG(ERROR) << StringPrintf("Ignoring splatted spectrum with negative "
                                   "luminance %f at (%f, %f)", v.y(), p.x, p.y);
        return;
    } else if (std::isinf(v.y())) {
        LOG(ERROR) << StringPrintf("Ignoring splatted spectrum with infinite "
                                   "luminance at (%f, %f)", p.x, p.y);
        return;
    }

    if (!InsideExclusive((Point2i)p, croppedPixelBounds)) return;
    if (v.y() > maxSampleLuminance)
        v *= maxSampleLuminance / v.y();
    Float xyz[3];
    v.ToXYZ(xyz);
    Pixel &pixel = GetPixel((Point2i)p);
    for (int i = 0; i < 3; ++i) pixel.splatXYZ[i].Add(xyz[i]);
}

void Film::UpdateDisplay(int x0, int y0, int x1, int y1, float splatScale) 
{
    if (!sdl_window) 
    {
		return;
	}
    // Compute window coordinates for pixels to update
    int xPixelStart =0;
    int xPixelCount = fullResolution.x;
    int yPixelStart = 0;
    int yPixelCount = fullResolution.y;
    x0 -= xPixelStart;
    x1 -= xPixelStart;
    y0 -= yPixelStart;
    y1 -= yPixelStart;
    x0 = Clamp(x0, 0, xPixelCount);
    x1 = Clamp(x1, 0, xPixelCount);
    y0 = Clamp(y0, 0, yPixelCount);
    y1 = Clamp(y1, 0, yPixelCount);
    u_int *pix = new u_int[(x1-x0)*(y1-y0)];
    u_int *pp = pix;
    for (int y = y0; y < y1; ++y) 
    {
        for (int x = x0; x < x1; ++x) 
        {
            // Compute weighted pixel value and update window pixel
            Point2i p(x,y);
            Pixel &pixel = GetPixel(p);
            Float rgb[3];
            XYZToRGB(pixel.xyz, rgb);
            Float weightSum = pixel.filterWeightSum;
            if (weightSum != 0.f) 
            {
                Float invWt = 1.f / weightSum;
                rgb[0] *= invWt;
                rgb[1] *= invWt;
                rgb[2] *= invWt;
            }

            Float splatRGB[3];
            Float splatXYZ[3] = {pixel.splatXYZ[0], pixel.splatXYZ[1], pixel.splatXYZ[2]};
            XYZToRGB(splatXYZ, splatRGB);
            rgb[0] += splatScale * splatRGB[0];
            rgb[1] += splatScale * splatRGB[1];
            rgb[2] += splatScale * splatRGB[2];
            
            *pp++ = SDL_MapRGB(sdl_surface->format,
                u_char(Clamp(powf(rgb[0], 1./1.8), 0.f, 1.f) * 255),
                u_char(Clamp(powf(rgb[1], 1./1.8), 0.f, 1.f) * 255),
                u_char(Clamp(powf(rgb[2], 1./1.8), 0.f, 1.f) * 255));
        }
    }
    // Acquire mutex lock, update window pixels, redraw
    //MutexLock lock(*sdlmutex);
    std::lock_guard<std::mutex> lock(sdl_mutex);
    if (SDL_LockSurface(sdl_surface) == -1) 
    { 
	}
    
    pp=pix;
    for (int y = y0; y < y1; ++y) 
    {
        for (int x = x0; x < x1; ++x) 
        {
            u_int *bufp = (u_int *)sdl_surface->pixels + y*sdl_surface->pitch/4 + x;
            *bufp = *pp++;
        }
    }
    SDL_UnlockSurface(sdl_surface);
    SDL_UpdateTexture(sdl_texture, NULL, sdl_surface->pixels, sdl_surface->pitch);
    SDL_RenderClear(sdl_renderer);
    SDL_RenderCopy(sdl_renderer, sdl_texture, NULL, NULL);
    SDL_RenderPresent(sdl_renderer);
    delete[] pix;
}


void Film::WriteImage(Float splatScale) {
    // Convert image to RGB and compute final pixel values
    LOG(INFO) <<
        "Converting image to RGB and computing final weighted pixel values";
    std::unique_ptr<Float[]> rgb(new Float[3 * croppedPixelBounds.Area()]);
    int offset = 0;
    for (Point2i p : croppedPixelBounds) {
        // Convert pixel XYZ color to RGB
        Pixel &pixel = GetPixel(p);
        XYZToRGB(pixel.xyz, &rgb[3 * offset]);

        // Normalize pixel with weight sum
        Float filterWeightSum = pixel.filterWeightSum;
        if (filterWeightSum != 0) {
            Float invWt = (Float)1 / filterWeightSum;
            rgb[3 * offset] = std::max((Float)0, rgb[3 * offset] * invWt);
            rgb[3 * offset + 1] =
                std::max((Float)0, rgb[3 * offset + 1] * invWt);
            rgb[3 * offset + 2] =
                std::max((Float)0, rgb[3 * offset + 2] * invWt);
        }

        // Add splat value at pixel
        Float splatRGB[3];
        Float splatXYZ[3] = {pixel.splatXYZ[0], pixel.splatXYZ[1],
                             pixel.splatXYZ[2]};
        XYZToRGB(splatXYZ, splatRGB);
        rgb[3 * offset] += splatScale * splatRGB[0];
        rgb[3 * offset + 1] += splatScale * splatRGB[1];
        rgb[3 * offset + 2] += splatScale * splatRGB[2];

        // Scale pixel value by _scale_
        rgb[3 * offset] *= scale;
        rgb[3 * offset + 1] *= scale;
        rgb[3 * offset + 2] *= scale;
        ++offset;
    }

    // Write RGB image
    LOG(INFO) << "Writing image " << filename << " with bounds " <<
        croppedPixelBounds;
    pbrt::WriteImage(filename, &rgb[0], croppedPixelBounds, fullResolution);
}

Film *CreateFilm(const ParamSet &params, std::unique_ptr<Filter> filter) {
    // Intentionally use FindOneString() rather than FindOneFilename() here
    // so that the rendered image is left in the working directory, rather
    // than the directory the scene file lives in.
    std::string filename = params.FindOneString("filename", "");
    if (PbrtOptions.imageFile != "") {
        if (filename != "") {
            Warning(
                "Output filename supplied on command line, \"%s\", ignored "
                "due to filename provided in scene description file, \"%s\".",
                PbrtOptions.imageFile.c_str(), filename.c_str());
        } else
            filename = PbrtOptions.imageFile;
    }
    if (filename == "") filename = "pbrt.exr";

    int xres = params.FindOneInt("xresolution", 1280);
    int yres = params.FindOneInt("yresolution", 720);
    if (PbrtOptions.quickRender) xres = std::max(1, xres / 4);
    if (PbrtOptions.quickRender) yres = std::max(1, yres / 4);
    Bounds2f crop(Point2f(0, 0), Point2f(1, 1));
    int cwi;
    const Float *cr = params.FindFloat("cropwindow", &cwi);
    if (cr && cwi == 4) {
        crop.pMin.x = Clamp(std::min(cr[0], cr[1]), 0.f, 1.f);
        crop.pMax.x = Clamp(std::max(cr[0], cr[1]), 0.f, 1.f);
        crop.pMin.y = Clamp(std::min(cr[2], cr[3]), 0.f, 1.f);
        crop.pMax.y = Clamp(std::max(cr[2], cr[3]), 0.f, 1.f);
    } else if (cr)
        Error("%d values supplied for \"cropwindow\". Expected 4.", cwi);

    Float scale = params.FindOneFloat("scale", 1.);
    Float diagonal = params.FindOneFloat("diagonal", 35.);
    Float maxSampleLuminance = params.FindOneFloat("maxsampleluminance",
                                                   Infinity);
    return new Film(Point2i(xres, yres), crop, std::move(filter), diagonal,
                    filename, scale, maxSampleLuminance);
}

}  // namespace pbrt
