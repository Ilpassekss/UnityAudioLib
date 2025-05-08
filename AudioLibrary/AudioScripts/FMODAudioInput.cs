using FMODUnity;
using UnityEngine;
using System;
using System.Runtime.InteropServices;
using System.Threading;
using System.Collections.Concurrent;
using SpectralAnalyze;
using System.IO;






public class FMODAudioInput : MonoBehaviour
{
    private FMOD.System fmodSystem;
    private FMOD.Sound recordingSound;
    private FMOD.Channel playbackChannel;
    private FMOD.RESULT result;
    public bool isRecording = false;
    [SerializeField] public bool playbackStarted = false;
    // for Android only one
    [SerializeField] WindowTypes windowType;
     [SerializeField]private int recordDeviceIndex = 0;
    // Samples buffer size
    [SerializeField] private uint playbackThreshold = 256;

    // FFT params
    [SerializeField] int sampleRate = 48000;
    [SerializeField] int fftSize = 8192; // ~11.71 (sampleRate/fftSize = freqStep)
    private FFTRealization fftProcessor;
    

    // Audio buffers
    private float[] audioBuffer;         
    private ComplexNumber[] fftOutput;   
    private short[] pcmBuffer;          
    private Thread audioProcessingThread;
    private bool processingThreadRunning = false;

    // data for shader
    [Header("Spectrum material")]
    [SerializeField] private Material spectrumMaterial;
    private float[] fftdB;
    private float[] fftdBCopy;// this info for shader
    private readonly object fftLock = new object();

    private ConcurrentQueue<string> logQueue = new ConcurrentQueue<string>();

    void Start()
    {
        fftProcessor = new FFTRealization(fftSize, fftSize, WindowTypes.eHammingWindow);
        audioBuffer = new float[fftSize];
        fftOutput = new ComplexNumber[fftSize];
        pcmBuffer = new short[fftSize];  

        //init buffers for shading
        fftdB = new float[fftSize/8];
        fftdBCopy = new float[fftSize/8];
        
        FMOD.Studio.System studioSystem = RuntimeManager.StudioSystem;
        studioSystem.getCoreSystem(out fmodSystem);

        int numDrivers = 0, numConnected = 0;
        result = fmodSystem.getRecordNumDrivers(out numDrivers, out numConnected);
        if (result != FMOD.RESULT.OK || numDrivers == 0)
        {
            Debug.LogError("Any audio input devices : " + result);
            return;
        }

        FMOD.CREATESOUNDEXINFO exinfo = new FMOD.CREATESOUNDEXINFO();
        exinfo.cbsize = Marshal.SizeOf(typeof(FMOD.CREATESOUNDEXINFO));
        exinfo.numchannels = 1;                // mono channel
        exinfo.defaultfrequency = sampleRate;  // 48000 Hz
        exinfo.format = FMOD.SOUND_FORMAT.PCM16; // 16-bit PCM
        exinfo.decodebuffersize = 256;         // low latency buffer
        exinfo.length = (uint)sampleRate * (uint)exinfo.numchannels * 2;  

       
        result = fmodSystem.createSound("", FMOD.MODE.LOOP_NORMAL | FMOD.MODE.OPENUSER, ref exinfo, out recordingSound);
        if (result != FMOD.RESULT.OK)
        {
            Debug.LogError("Error: " + result);
            return;
        }

        
        result = fmodSystem.recordStart(recordDeviceIndex, recordingSound, true);
        if (result != FMOD.RESULT.OK)
        {
            Debug.LogError("Recording start error: " + result);
            return;
        }
        isRecording = true;
        
        Debug.Log("Micro recording started, waiting for data...");

        processingThreadRunning = true;
        audioProcessingThread = new Thread(AudioProcessingThread);
        audioProcessingThread.Start();
    }

    void Update()
    {
        if (!isRecording)
            return;

        
        if (!playbackStarted)
        {
            uint recPos = 0;
            result = fmodSystem.getRecordPosition(recordDeviceIndex, out recPos);
            if (result != FMOD.RESULT.OK)
                return;

            if (recPos > playbackThreshold)
            {
                //playing sound back in speakers 
                // result = fmodSystem.playSound(recordingSound, default(FMOD.ChannelGroup), false, out playbackChannel);
                // if (result != FMOD.RESULT.OK)
                // {
                //     Debug.LogError("Ошибка при запуске воспроизведения: " + result);
                //     return;
                // }
                playbackStarted = true;
                Debug.Log("Воспроизведение запущено.");
            }
            
            lock(fftLock)
            {
                
                spectrumMaterial.SetFloatArray("", fftdBCopy);
            }

        }

        while (logQueue.TryDequeue(out string message))
        {
            Debug.Log(message);
            
        }
    }

    // Audio 
    private void AudioProcessingThread()
    {
        logQueue.Enqueue("Thread started");

        while (processingThreadRunning)
        {
            if (playbackStarted)
            {
                if (GetCurrentAudioData())
                {
                    fftProcessor.DFT(audioBuffer, fftOutput, 1);

                    
                    int maxBin = 0;
                    float maxMagnitude = 0f;
                    
                    lock(fftLock)
                    {
                        float magnitudesSumForNormalization = 0f;

                        for (int i = 1; i < fftSize / 2; i++)
                        {

                            float magnitude = (float)fftOutput[i].Magnitude;
                            float safeMagnitude = Mathf.Max(magnitude, 1e-7f);

                            if((float)i % 8f != 0f)
                            {
                                magnitudesSumForNormalization += safeMagnitude;
                            }
                            else
                            {
                                fftdB[i-7] = 20f * Mathf.Log10(magnitudesSumForNormalization / 8f ); 

                                magnitudesSumForNormalization = 0f;
                            }

                            // fftdB[i] = 20f * Mathf.Log10(safeMagnitude);

                            if (magnitude > maxMagnitude)
                            {
                                maxMagnitude = magnitude;
                                maxBin = i;
                            }
                        }
                        Array.Copy(fftdB, fftdBCopy, fftdB.Length);
                    }

                    // check interpolation opportunity
                    if (maxBin > 0 && maxBin < (fftSize / 2) - 1)
                    {
                        // magnitudes for interpolation
                        float magL = (float)fftOutput[maxBin - 1].Magnitude;
                        float magC = (float)fftOutput[maxBin].Magnitude;
                        float magR = (float)fftOutput[maxBin + 1].Magnitude;

                        // safe interpolation
                        float denom = (magL - 2f * magC + magR);
                        float bitCorrection = 0f;
                        if (Mathf.Abs(denom) > 1e-6f)
                        {
                            bitCorrection = 0.5f * (magL - magR) / denom;
                        }

                        float interpolatedBin = maxBin + bitCorrection;
                        float freq = interpolatedBin * sampleRate / fftSize;

                        // уровень сигнала пересчитываем по центру
                        float interpolatedMagnitude = magC; // или можно взять с учетом интерполяции — но можно и просто magC
                        float safeMagnitude = Mathf.Max(interpolatedMagnitude, 1e-7f);
                        float dB = 20f * Mathf.Log10(safeMagnitude);

                        if (dB > -40f)
                        {
                            string logLine = $"Peak frequency: {freq:F2} Hz, Level: {dB:F2} dB";
                            File.AppendAllText("Assets/Scripts/AudioScripts/audiologs.txt", logLine + Environment.NewLine);
                        }
                    }
                    else
                    {
                        // Без интерполяции (края спектра)
                        float freq = maxBin * sampleRate / fftSize;
                        float safeMagnitude = Mathf.Max(maxMagnitude, 1e-7f);
                        float dB = 20f * Mathf.Log10(safeMagnitude);

                        if (dB > -40f)
                        {
                            string logLine = $"Peak frequency: {freq:F2} Hz, Level: {dB:F2} dB (no interpolation)";
                            File.AppendAllText("Assets/Scripts/AudioScripts/audiologs.txt", logLine + Environment.NewLine);
                        }
                    }


                }
            }
            // 
            Thread.Sleep(50);
        }
    }

    bool GetCurrentAudioData()
    {
        uint recPos;
        result = fmodSystem.getRecordPosition(recordDeviceIndex, out recPos);
        if (result != FMOD.RESULT.OK)
        {
            logQueue.Enqueue("Ошибка при получении позиции записи.");
            return false;
        }

        IntPtr ptr1, ptr2;
        uint len1, len2;
        uint soundLength;

        recordingSound.getLength(out soundLength, FMOD.TIMEUNIT.PCM);

        
        // fmod can work like circled buffer 
        uint startPos = (recPos >= fftSize) ? (recPos - (uint)fftSize) : 0;

        // buffer locking for reading
        result = recordingSound.@lock(startPos, (uint)fftSize, out ptr1, out ptr2, out len1, out len2);

        if (result != FMOD.RESULT.OK)
        {
            logQueue.Enqueue("Sound buffer locking error.");
            return false;
        }

        // copy data into audio buffer from block 1 & 2 in PCM format into pcmBuffer
        if (len1 > 0)
            Marshal.Copy(ptr1, pcmBuffer, 0, (int)len1 / 2);
        if (len2 > 0)
            Marshal.Copy(ptr2, pcmBuffer, (int)len1 / 2, (int)len2 / 2);

        // buffer unlocking
        recordingSound.unlock(ptr1, ptr2, len1, len2);

        // samples counting 
        int totalSamples = (int)((len1 + len2) / 2);
        // copy data in in buffer 
        for (int i = 0; i < fftSize; i++)
        {
            // normlizing [-32768, 32767] diapazone in [-1.0, 1.0] format
            audioBuffer[i] = (i < totalSamples) ? pcmBuffer[i] / 32768.0f : 0.0f;
        }

        return true;
    }

    void OnDestroy()
    {
        processingThreadRunning = false;
        if (audioProcessingThread != null && audioProcessingThread.IsAlive)
        {
            audioProcessingThread.Join();
        }

        if (isRecording)
        {
            fmodSystem.recordStop(recordDeviceIndex);
            isRecording = false;
        }
        if (recordingSound.handle != IntPtr.Zero)
        {
            recordingSound.release();
        }
    }
}

