using System;
using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using UnityEngine;


namespace SpectralAnalyze
{

    public enum WindowTypes
    {
        eRectangular = 0,
        eHammingWindow,
        eHanningWindow,
        eBlackmanWindow,
        eBlackmanHarrisWindow,
        
        // Additional common windows
        eBartlettWindow,        // Triangular window
        eFlatTopWindow,
        eGaussianWindow,
        eKaiserWindow,
        eTukeyWindow,           // Tapered cosine
        
        // Advanced window functions
        eDolphChebyshevWindow,
        eNuttallWindow,
        eParzenWindow,
        eWelchWindow,
        eBohmanWindow,
        
        // Specialized windows
        ePoissonWindow,
        eLanczosWindow,
        eExponentialWindow,
        eUltrasphericalWindow,
        
        // Blackman variations
        eBlackmanNuttallWindow,
        eBlackmanExactWindow,
        
        // Other specialized windows
        eCosineWindow,          // Single raised cosine lobe
        eRiessBesselWindow,
        eBartlettHannWindow,
        ePlanchetWindow
    }

    public class FFTRealization
    {

        // arr of singnals in complex numbers form
        ComplexNumber[] mBuffer;
        // pre. calc. sin and cos values
        // make calc. faster 
        ComplexNumber[] mSinCosTable;
        // frame length
        int mFrameLength;
        // fft transf. length. Size (256, 521, etc)
        int mFftLength;
        // calc like log2(mFftLength)
        // mean count of iterations in FFT alg. 
        int mFftPower;

        // coof of window function
        // save from specktrum leakage
        float[] m_window;
        WindowTypes windowType;

        int[] mBitReversedIndexes;

        public FFTRealization(int fftLength)
        {
            this.Init(fftLength, fftLength, WindowTypes.eHammingWindow);
        }

        public FFTRealization(int fftLength, int frameLength)
        {
            this.Init(fftLength, frameLength, WindowTypes.eHammingWindow);
        }

        public FFTRealization(int fftLength, int frameLength, WindowTypes windowType)
        {
            this.Init(fftLength, frameLength, windowType);
        }

        public void Init(int mFftLength, int frameLength, WindowTypes windowType)
        {

            double dPower = Math.Log(Convert.ToDouble(mFftLength) / Math.Log(2.0));

            this.windowType = windowType;
            this.mFrameLength = frameLength;
            this.mFftLength = mFftLength;
            this.mFftPower = (int)dPower;
            this.mSinCosTable = new ComplexNumber[mFftLength];
            this.mBitReversedIndexes = new int[mFftLength];
            this.m_window = new float[frameLength];
            this.mBuffer = new ComplexNumber[frameLength];


            // precalc sin cos lookup table
            double ang = 2.0 * Math.PI / Convert.ToDouble(mFftLength);
            for(int i = 0; i < mFftLength; i ++ )
            {
                mSinCosTable[i].x = (float)Math.Cos(ang * Convert.ToDouble(i));
                mSinCosTable[i].y = (float)Math.Sin(ang * Convert.ToDouble(i));
            }

            // create bit reversed indexes
            for (int k = 0; k < mFftLength; k++)
            {
                int bitL = mFftLength >> 1;
                int C = 0;
                int tmp = k;

                for (int i = 0; i < mFftPower; i++)
                {
                    if ((tmp & 1) > 0)
                        C += bitL;
                    tmp >>= 1;
                    bitL >>= 1;
                }

                mBitReversedIndexes[k] = C;
            }   

            // create Hamming window
            float a0 = 0.0f;
            float a1 = 0.0f;
            float a2 = 0.0f;
            float a3 = 0.0f;
            float angle1 = (float)(2.0 * Math.PI / Convert.ToDouble(frameLength - 1));
            float angle2 = (float)(4.0 * Math.PI / Convert.ToDouble(frameLength - 1));
            float angle3 = (float)(6.0 * Math.PI / Convert.ToDouble(frameLength - 1));  

            switch(windowType)
            {
                case WindowTypes.eRectangular:
                    a0 = 1.0f;
                    break;

                case WindowTypes.eHammingWindow:
                    a0 = 0.54f;
                    a1 = 0.46f;
                    break;

                case WindowTypes.eHanningWindow:
                    a0 = 0.5f;
                    a1 = 0.5f;
                    break;

                case WindowTypes.eBlackmanWindow:
                    a0 = 0.42659f;
                    a1 = 0.49656f;
                    a2 = 0.076849f;
                    break;

                case WindowTypes.eBlackmanHarrisWindow:
                    a0 = 0.35875f;
                    a1 = 0.48829f;
                    a2 = 0.14128f;
                    a3 = 0.01168f;
                    break;                
            }   

            for (int i = 0; i < frameLength; i++)
            {
                double result = a0 - a1 * Math.Cos(angle1 * Convert.ToDouble(i)) + a2 * Math.Cos(angle2 * Convert.ToDouble(i)) - a3 * Math.Cos(angle3 * Convert.ToDouble(i));
                m_window[i] = (float)result;
            }    

        }

        public float[,] GetWindowResponse(int frameSize, int fftSize, int oversampleRatio, int interpolationSize, WindowTypes windowType)
        {
            
            fftSize *= oversampleRatio;

            int offset = -( interpolationSize / 2 ) * oversampleRatio;
            int mask = fftSize - 1;
            FFTRealization ffTRealization = new FFTRealization(fftSize, frameSize, windowType);
            ComplexNumber[] fftBuffer = new ComplexNumber[fftSize];
            float[,] outputBuffer = new float[oversampleRatio, interpolationSize];

            ffTRealization.FFT(fftBuffer);

            // DSPHelpers.NormalizePeak(fftBuffer);
            for(int i = 0, index = offset; i < oversampleRatio; i++, index --)
            {
                for(int j = 0 , index2 = index; j< interpolationSize; j++, index2 += oversampleRatio)
                {
                    outputBuffer[i, j] = (float)fftBuffer[index2 & mask].Magnitude;
                }
            }


            return outputBuffer;
        }
        
        public void IDFT(ComplexNumber[] input, float[] output)
        {
            float gain = 1.0f;
            float maxVal = 0.0f;
            bool bNormalize = false;

            for (int k = 0; k < mFftLength; k += 1)
            {
                float sum;
                int freq = 0;
                int mask = mFftLength - 1;

                sum = 0.0f;

                for (int n = 0; n < mFftLength; n++)
                {
                    sum += input[n].x * mSinCosTable[freq].x - input[n].y * mSinCosTable[freq].y;
                    freq = (freq + k);

                    if (freq >= mFftLength)
                        freq -= mFftLength;
                }

                output[k] = sum * gain;
                if (Math.Abs(output[k]) > maxVal)
                    maxVal = Math.Abs(output[k]);
            }

            if (bNormalize == true)
            {
                gain = 1.0F / maxVal;
                for (int i = 0; i < mFftLength; i++)
                {
                    output[i] *= gain;
                }
            }
        }

        public void DFT(float[] input, ComplexNumber[] output, int interval)
        {
            WindowInput(input, mBuffer);
            int lenOver2 = mFftLength >> 1;

            for (int k = 0; k < mFftLength; k += interval)
            {
                
                int freq = 0;
                int mask = mFftLength - 1;

                ComplexNumber sum = new ComplexNumber(0.0f, 0.0f);

                for (int n = 0; n < mFftLength; n++)
                {
                    sum.x += input[n] * mSinCosTable[freq].x;
                    sum.y -= input[n] * mSinCosTable[freq].y;

                    freq = freq + k;

                    if (freq >= mFftLength)
                        freq -= mFftLength;

                }

                output[k] = sum;
            }
            
        }

        public void FFT(float[] input, ComplexNumber[] output)
        {
            WindowInput(input, mBuffer);
            DecimationInFrequencyFFT(mBuffer, output);
        }

        public void FFT(float[] input, int offset, ComplexNumber[] output)
        {
            WindowInput(input, offset, mBuffer);
            DecimationInFrequencyFFT(mBuffer, output);
        }

        public void FFT(ComplexNumber[] output)
        {
            WindowInput(null, mBuffer);
            DecimationInFrequencyFFT(mBuffer, output);
        }

        public void FFT(ComplexNumber[] input, ComplexNumber[] output)
        {
            DecimationInFrequencyFFT(input, output);
        }

        private void WindowInput(float[] input, ComplexNumber[] output)
        {
            // this function multiplies the time signals with the window signal
            // padds zeros to the end of the sequences if necessary
            float scaler = 1.0f / Convert.ToSingle(mFrameLength);

            if (input != null)
            {
                for (int i = 0; i < mFrameLength; i++)
                {
                    output[i].x = scaler * input[i] * m_window[i];
                    output[i].y = 0.0F;
                }
            }
            else
            {
                for (int i = 0; i < mFrameLength; i++)
                {
                    output[i].x = scaler * m_window[i];
                    output[i].y = 0.0F;
                }
            }

            // zero pad the rest
            for (int i = mFrameLength; i < mFrameLength; i++)
            {
                output[i].x = 0.0f;
                output[i].y = 0.0f;
            }
        }

        private void WindowInput(float[] input, int offset, ComplexNumber[] output)
        {
            // this function multiplies the time signals with the window signal
            // padds zeros to the end of the sequences if necessary
            float scaler = 1.0f / Convert.ToSingle(mFrameLength);
            int inputLength = input.Length;

            for (int i = 0, j = offset; i < mFrameLength; i++, j++)
            {
                float data = (j < inputLength) ? input[j] : 0.0f;
                output[i].x = scaler * data * m_window[i];
                output[i].y = 0.0F;
            }

            // zero pad the rest
            for (int i = mFrameLength; i < mFrameLength; i++)
            {
                output[i].x = 0.0f;
                output[i].y = 0.0f;
            }
        }
    

        private void DecimationInFrequencyFFT(ComplexNumber[] Input, ComplexNumber[] output)
        {
            // this fft algorithm uses decimation in frequency
            // the input may come in bit reversed index form
            // n is the FFT length
            // m is the FFT power -> log2(n);
            int n = mFftLength;
            int m = mFftPower - 1;
            int indexIncrement = n;
            int innerLoops = 1;
            int innerLoopLength = n >> 1;

            for (int iter = 0; iter < m; iter++)
            {
                for (int index = 0, loops = 0; loops < innerLoops; loops++, index += indexIncrement)
                {
                    for (int i = index, j = 0, angle = 0, angleIncrement = innerLoops; j < innerLoopLength; j++, i++, angle += angleIncrement)
                    {
                        int k = i + innerLoopLength;
                        ComplexNumber T = Input[k];
                        float x = Input[i].x - T.x;
                        float y = Input[i].y - T.y;
                        Input[i].x += T.x;
                        Input[i].y += T.y;
                        Input[k].x = x * mSinCosTable[angle].x + y * mSinCosTable[angle].y;
                        Input[k].y = y * mSinCosTable[angle].x - x * mSinCosTable[angle].y;
                    }
                }

                indexIncrement >>= 1;
                innerLoops <<= 1;
                innerLoopLength >>= 1;
            }

            for (int i = 0; i < n; i += 2)
            {
                int j = i + 1;
                int reverseI = mBitReversedIndexes[i];
                int reverseJ = mBitReversedIndexes[j];
                ComplexNumber InputI = Input[i];
                ComplexNumber InputJ = Input[j];
                output[reverseI].x = Input[i].x + Input[j].x;
                output[reverseI].y = Input[i].y + Input[j].y;
                output[reverseJ].x = Input[i].x - Input[j].x;
                output[reverseJ].y = Input[i].y - Input[j].y;
            }

            for (int i = 0; i < n; i++)
            {
                if (output[i].Magnitude < 0.000001)
                {
                    output[i].x = 0.0f;
                    output[i].y = 0.0f;
                }
            }
        }

        private void DecimationInFrequencyIFFT(ComplexNumber[] Input, float[] output)
        {
            // this fft algorithm uses decimation in frequency
            // the input may come in bit reversed index form
            // n is the FFT length
            // m is the FFT power -> log2(n);
            int n = mFftLength;
            int m = mFftPower - 1;
            int indexIncrement = n;
            int innerLoops = 1;
            int innerLoopLength = n >> 1;

            for (int iter = 0; iter < m; iter++)
            {
                for (int index = 0, loops = 0; loops < innerLoops; loops++, index += indexIncrement)
                {
                    for (int i = index, j = 0, angle = 0, angleIncrement = innerLoops; j < innerLoopLength; j++, i++, angle += angleIncrement)
                    {
                        int k = i + innerLoopLength;
                        ComplexNumber T = Input[k];
                        float x = Input[i].x - T.x;
                        float y = Input[i].y - T.y;
                        Input[i].x += T.x;
                        Input[i].y += T.y;
                        Input[k].x = x * mSinCosTable[angle].x - y * mSinCosTable[angle].y;
                        Input[k].y = y * mSinCosTable[angle].x + x * mSinCosTable[angle].y;
                    }
                }

                indexIncrement >>= 1;
                innerLoops <<= 1;
                innerLoopLength >>= 1;
            }

            for (int i = 0; i < n; i += 2)
            {
                int j = i + 1;
                int reverseI = mBitReversedIndexes[i];
                int reverseJ = mBitReversedIndexes[j];
                output[reverseI] = Input[i].x + Input[j].x;
                output[reverseJ] = Input[i].x - Input[j].x;
            }
        }
        private void DecimationInFrequencyIFFT(ComplexNumber[] Input, float[] output, int offset)
        {
            // this fft algorithm uses decimation in frequency
            // the input may come in bit reversed index form
            // n is the FFT length
            // m is the FFT power -> log2(n);
            int n = mFftLength;
            int m = mFftPower - 1;
            int indexIncrement = n;
            int innerLoops = 1;
            int innerLoopLength = n >> 1;

            for (int iter = 0; iter < m; iter++)
            {
                for (int index = 0, loops = 0; loops < innerLoops; loops++, index += indexIncrement)
                {
                    for (int i = index, j = 0, angle = 0, angleIncrement = innerLoops; j < innerLoopLength; j++, i++, angle += angleIncrement)
                    {
                        int k = i + innerLoopLength;
                        ComplexNumber T = Input[k];
                        float x = Input[i].x - T.x;
                        float y = Input[i].y - T.y;
                        Input[i].x += T.x;
                        Input[i].y += T.y;
                        Input[k].x = x * mSinCosTable[angle].x - y * mSinCosTable[angle].y;
                        Input[k].y = y * mSinCosTable[angle].x + x * mSinCosTable[angle].y;
                    }
                }

                indexIncrement >>= 1;
                innerLoops <<= 1;
                innerLoopLength >>= 1;
            }

            for (int i = 0; i < n; i += 2)
            {
                int j = i + 1;
                int reverseI = offset + mBitReversedIndexes[i];
                int reverseJ = offset + mBitReversedIndexes[j];
                output[reverseI] = Input[i].x + Input[j].x;
                output[reverseJ] = Input[i].x - Input[j].x;
            }
        }    
    
    }
}
