#include <SDL2/SDL.h>
#include <math.h>
#include <stdio.h>
#include <climits>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cstring>
#include "FastNoiseLite.h"
#include <chrono>
#include <mutex>
#include <thread>

#define SCREEN_WIDTH 500
#define SCREEN_HEIGHT 300   
#define NOISE_FREQ 1
#define RED_FREQ 10
#define GREEN_FREQ 50
#define BLUE_FREQ 1
#define OCTAVES 16
 
static std::mutex mtx_cout;

struct acout
{
        std::unique_lock<std::mutex> lk;
        acout()
            :
              lk(std::unique_lock<std::mutex>(mtx_cout))
        {

        }

        template<typename T>
        acout& operator<<(const T& _t)
        {
            std::cout << _t;
            return *this;
        }

        acout& operator<<(std::ostream& (*fp)(std::ostream&))
        {
            std::cout << fp;
            return *this;
        }
};

void updateNoiseField(uint* i, uint* j, float *xred_off, float *xgreen_off, float *yred_off, float *ygreen_off, std::vector<unsigned char>* pixels, std::vector<unsigned char>* pixels2, FastNoiseLite* noise, size_t threads) {
                for (int repeat = 0; repeat < (SCREEN_WIDTH-1)*((SCREEN_HEIGHT-2)/(int)(threads-1)); repeat++) {
                *i += 1;
                if (*i < SCREEN_WIDTH-1) {
                    if (*j > 0) {
                        float xfactor = (float)*i/(float)SCREEN_WIDTH;
                        float yfactor = (float)*j/(float)SCREEN_HEIGHT;

                        const unsigned int offset = ( SCREEN_WIDTH * 4 * *j ) + *i * 4;

                        const unsigned int leftneighbour = ( SCREEN_WIDTH * 4 * (*j+1) ) + (*i-1) * 4; 
                        const unsigned int rightneighbour = ( SCREEN_WIDTH * 4 * (*j+1) ) + (*i+1) * 4; 
                        const unsigned int topneighbour = ( SCREEN_WIDTH * 4 * (*j-1) ) + *i * 4; 
                        const unsigned int bottomneighbour = ( SCREEN_WIDTH * 4 * (*j+1) ) + *i * 4;
                        uint8_t tmp_noiseR = (uint8_t)(floor(255 * noise->GetNoise((xfactor + *xred_off)*RED_FREQ, (yfactor + *yred_off)*RED_FREQ)));
                        uint8_t tmp_noiseG = (uint8_t)(floor(128 * noise->GetNoise((xfactor + *xgreen_off)*GREEN_FREQ, (yfactor + *ygreen_off)*GREEN_FREQ)));

                        tmp_noiseR /= floor(100+(rand()%100)+(*j/SCREEN_HEIGHT));
                        tmp_noiseG /= floor(25+rand()%230);
                        if (*i > 1 && *i < SCREEN_WIDTH-2 && *j > 1 && *j < SCREEN_HEIGHT-2) {
                                int tmp_pixelB = 0;
                                int tmp_pixelR = ((*pixels2)[leftneighbour+2] + (*pixels2)[rightneighbour+2] + (*pixels2)[topneighbour+2] + (*pixels2)[bottomneighbour+2])/4 - tmp_noiseR;
                                int tmp_pixelG = ((*pixels2)[leftneighbour+1] + (*pixels2)[rightneighbour+1] + (*pixels2)[topneighbour+1] + (*pixels2)[bottomneighbour+1])/4 - tmp_noiseG;
                                if (tmp_pixelR <= 0)
                                    tmp_pixelR = 0;
                                if (tmp_pixelR >= 255)
                                    tmp_pixelR = 255;
                                if (tmp_pixelG <= 0)
                                    tmp_pixelG = 0;
                                if (tmp_pixelG >= 255)
                                    tmp_pixelG = 255;
                                (*pixels)[ offset + 0 ] = (uint8_t)tmp_pixelB;
                                (*pixels)[ offset + 1 ] = (uint8_t)tmp_pixelG;
                                (*pixels)[ offset + 2 ] = (uint8_t)tmp_pixelR;
                                (*pixels)[ offset + 3 ] = SDL_ALPHA_OPAQUE;
                        }
                    }
                    else 
                    {
                        *j = SCREEN_HEIGHT-3;
                        //*yblue_off += 0.02;
                        *ygreen_off += 0.15;
                        *yred_off += 0.20;
                    }
                
                }
                    else {
                        *i = 0;
                        *j -= 1;
                    }
                }
            
}

int main( int argc, char** argv )
{
    SDL_Init( SDL_INIT_EVERYTHING );

    SDL_Window* window = SDL_CreateWindow
        (
        "Noise",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
        SCREEN_WIDTH, SCREEN_HEIGHT,
        SDL_WINDOW_FULLSCREEN_DESKTOP
        );

    SDL_Renderer* renderer = SDL_CreateRenderer
        (
        window,
        -1,
        SDL_RENDERER_ACCELERATED
        );

    SDL_RendererInfo info;
    SDL_GetRendererInfo( renderer, &info );
    std::cout << "Renderer name: " << info.name << std::endl;
    std::cout << "Texture formats: " << std::endl;
    for( Uint32 i = 0; i < info.num_texture_formats; i++ )
    {
        std::cout << SDL_GetPixelFormatName( info.texture_formats[i] ) << std::endl;
    }

    const unsigned int texWidth = SCREEN_WIDTH;
    const unsigned int texHeight = SCREEN_HEIGHT;
    SDL_Texture* texture = SDL_CreateTexture
        (
        renderer,
        SDL_PIXELFORMAT_ARGB8888,
        SDL_TEXTUREACCESS_STREAMING,
        texWidth, texHeight
        );

    std::vector< unsigned char > pixels( texWidth * texHeight * 4, 0 );
    std::vector< unsigned char > pixels2( texWidth * texHeight * 4, 0 );

    SDL_Event event;
    bool running = true;
    bool useLocktexture = false;
    
	FastNoiseLite noise(1337);
	noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
	noise.SetFractalType(FastNoiseLite::FractalType_DomainWarpIndependent);
	noise.SetFractalOctaves(OCTAVES);
    noise.SetFrequency(NOISE_FREQ);

    unsigned int frames = 0;
	uint i = 0, j = 0;
    float xred_off = 0, xgreen_off = 100, prev_xred_off = 256, prev_xblue_off = 256, prev_xgreen_off = 256, yred_off = 0, yblue_off = 0, ygreen_off = 0, prev_yred_off = 256, prev_yblue_off = 256, prev_ygreen_off = 256;

    size_t thread_counter = 0;

    /*for (i = 0; i < SCREEN_WIDTH; i++) {
        for (j = 0; j < SCREEN_HEIGHT; j++) {

                const unsigned int offset = ( texWidth * 4 * j ) + i * 4;//= i * 4 * SCREEN_HEIGHT + j * 4; 
                pixels[ offset + 0 ] = 0;//(uint8_t)(floor(255 * noise. GetNoise((x)* BLUE_FREQ, (y)*BLUE_FREQ)));        // b
                pixels[ offset + 1 ] = 0;//(uint8_t)(floor(255 * noise.GetNoise((xfactor + xgreen_off)*GREEN_FREQ , (yfactor + xgreen_off)*GREEN_FREQ)));        // g
                pixels[ offset + 2 ] = 0;//(uint8_t)(floor(255 * noise.GetNoise((xfactor + xred_off)*RED_FREQ, (yfactor + xred_off)*RED_FREQ)));        // r
                pixels[ offset + 3 ] = SDL_ALPHA_OPAQUE;    // a
                pixels2[ offset + 0 ] = 0;//(uint8_t)(floor(255 * noise. GetNoise((x)* BLUE_FREQ, (y)*BLUE_FREQ)));        // b
                pixels2[ offset + 1 ] = 0;//(uint8_t)(floor(255 * noise.GetNoise((xfactor + xgreen_off)*GREEN_FREQ , (yfactor + xgreen_off)*GREEN_FREQ)));        // g
                pixels2[ offset + 2 ] = 0;//(uint8_t)(floor(255 * noise.GetNoise((xfactor + xred_off)*RED_FREQ, (yfactor + xred_off)*RED_FREQ)));        // r
                pixels2[ offset + 3 ] = SDL_ALPHA_OPAQUE;    // a
            }
        }*/
    j = SCREEN_HEIGHT-2;
    for (i = 0; i < SCREEN_WIDTH; i++) {

          const unsigned int offset = ( texWidth * 4 * j ) + i * 4;//= i * 4 * SCREEN_HEIGHT + j * 4; 
                pixels[ offset + 0 ] = 0;//(uint8_t)(floor(255 * noise. GetNoise((x)* BLUE_FREQ, (y)*BLUE_FREQ)));        // b
                pixels[ offset + 1 ] = 128;//(uint8_t)(floor(255 * noise.GetNoise((xfactor + xgreen_off)*GREEN_FREQ , (yfactor + xgreen_off)*GREEN_FREQ)));        // g
                pixels[ offset + 2 ] = 255;//(uint8_t)(floor(255 * noise.GetNoise((xfactor + xred_off)*RED_FREQ, (yfactor + xred_off)*RED_FREQ)));        // r
                pixels[ offset + 3 ] = SDL_ALPHA_OPAQUE;    // a
    } 
    j = SCREEN_HEIGHT-1;
    for (i = 0; i < SCREEN_WIDTH; i++) {

          const unsigned int offset = ( texWidth * 4 * j ) + i * 4;//= i * 4 * SCREEN_HEIGHT + j * 4; 
                pixels[ offset + 0 ] = 0;//(uint8_t)(floor(255 * noise. GetNoise((x)* BLUE_FREQ, (y)*BLUE_FREQ)));        // b
                pixels[ offset + 1 ] = 128;//(uint8_t)(floor(255 * noise.GetNoise((xfactor + xgreen_off)*GREEN_FREQ , (yfactor + xgreen_off)*GREEN_FREQ)));        // g
                pixels[ offset + 2 ] = 255;//(uint8_t)(floor(255 * noise.GetNoise((xfactor + xred_off)*RED_FREQ, (yfactor + xred_off)*RED_FREQ)));        // r
                pixels[ offset + 3 ] = SDL_ALPHA_OPAQUE;    // a
    }
    j = SCREEN_HEIGHT;
    for (i = 0; i < SCREEN_WIDTH; i++) {

          const unsigned int offset = ( texWidth * 4 * j ) + i * 4;//= i * 4 * SCREEN_HEIGHT + j * 4; 
                pixels[ offset + 0 ] = 0;//(uint8_t)(floor(255 * noise. GetNoise((x)* BLUE_FREQ, (y)*BLUE_FREQ)));        // b
                pixels[ offset + 1 ] = 128;//(uint8_t)(floor(255 * noise.GetNoise((xfactor + xgreen_off)*GREEN_FREQ , (yfactor + xgreen_off)*GREEN_FREQ)));        // g
                pixels[ offset + 2 ] = 255;//(uint8_t)(floor(255 * noise.GetNoise((xfactor + xred_off)*RED_FREQ, (yfactor + xred_off)*RED_FREQ)));        // r
                pixels[ offset + 3 ] = SDL_ALPHA_OPAQUE;    // a
    }
    //j = SCREEN_HEIGHT-1;
     

    Uint64 start = SDL_GetPerformanceCounter();

    while( running )
    {

        SDL_SetRenderDrawColor( renderer, 0, 0, 0, SDL_ALPHA_OPAQUE );
        SDL_RenderClear( renderer );

        while( SDL_PollEvent( &event ) )
        {
            if( ( SDL_QUIT == event.type ) ||
                ( SDL_KEYDOWN == event.type && SDL_SCANCODE_ESCAPE == event.key.keysym.scancode ) )
            {
                running = false;
                break;
            }
            if( SDL_KEYDOWN == event.type && SDL_SCANCODE_L == event.key.keysym.scancode )
            {
                useLocktexture = !useLocktexture;
                std::cout << "Using " << ( useLocktexture ? "SDL_LockTexture() + memcpy()" : "SDL_UpdateTexture()" ) << std::endl;
            }
        }
        
        //std::vector<std::thread> workers_cout;
        std::vector<std::thread> workers_acout;

        size_t worker(0);
        size_t threads(10);

        //std::cout << "\nWith acout():" << std::endl;
        //size_t tmp = thread_counter;
        for (int tmpint = 0; tmpint < 2; tmpint++) {
            if (thread_counter < threads) {
                //size_t tmp = thread_counter;
                //for (thread_counter; thread_counter < tmp+3; thread_counter++) { 
                i = 0;
                j = (uint)(SCREEN_HEIGHT/(uint)(threads)*(uint)thread_counter);
                //std::cout << j << std::endl;
                    workers_acout.emplace_back([&]
                        {
                            updateNoiseField(&i, &j, &xred_off, &xgreen_off, &yred_off, &ygreen_off, &pixels, &pixels2, &noise, threads);
                        }); 
                //worker++;
                //}
                thread_counter += 1;
                //++worker;
            } else {
                thread_counter = 0;
                //acout() << "\tThis is worker " << ++worker << " in thread " << std::this_thread::get_id() << std::endl;
            }
        }
        /*for (auto& w : workers_acout) {
            w.join();
        }*/
        //if (worker >= threads) {
            for (auto& w : workers_acout) {
                w.join();
            }
            worker = 0;
        //}

        
        

        if( useLocktexture )
        {
            unsigned char* lockedPixels = nullptr;
            int pitch = 0;
            SDL_LockTexture
                (
                texture,
                NULL,
                reinterpret_cast< void** >( &lockedPixels ),
                &pitch
                );
            std::memcpy(lockedPixels, pixels.data(), pixels.size());
            SDL_UnlockTexture( texture );
        }
        else
        {
            SDL_UpdateTexture
                (
                texture,
                NULL,
                pixels.data(),
                texWidth * 4
                );
        }

        pixels2 = pixels;
        SDL_RenderCopy( renderer, texture, NULL, NULL );
        SDL_RenderPresent( renderer );
        
        frames++;
        const Uint64 end = SDL_GetPerformanceCounter();
        const static Uint64 freq = SDL_GetPerformanceFrequency();
        const double seconds = ( end - start ) / static_cast< double >( freq );
        if( seconds > 2.0 )
        {
            std::cout
                << frames << " frames in "
                << std::setprecision(1) << std::fixed << seconds << " seconds = "
                << std::setprecision(1) << std::fixed << frames / seconds << " FPS ("
                << std::setprecision(3) << std::fixed << ( seconds * 1000.0 ) / frames << " ms/frame)"
                << std::endl;
            start = end;
            frames = 0;
        }
    }

    SDL_DestroyRenderer( renderer );
    SDL_DestroyWindow( window );
    SDL_Quit();

    return 0;
}
