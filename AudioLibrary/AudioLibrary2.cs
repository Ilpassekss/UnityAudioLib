using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Linq.Expressions;
using System.Numerics;
using UnityEngine;
using UnityEngine.UIElements;



namespace AudioLibrary
{
    public class SaveAudioClip
    {
        private const int HEADER_SIZE = 44;

        public static bool saveAudioClip(string filename, AudioClip audioClip)
        {

            if (!filename.ToLower().EndsWith(".wav"))
            {
                filename += ".wav";
            }

            var filepath = Path.Combine(Application.persistentDataPath, filename);// change filepath to android filepath

            Debug.Log(filepath);

            Directory.CreateDirectory(Path.GetDirectoryName(filepath));

            using (var fileStream = CreateEmpty(filepath))
            {

                //var trimmedClip = TrimSilence(audioClip, 0f);

                ConvertAndWrite(fileStream, audioClip);

                WriteHeader(fileStream, audioClip);
            }

            return true;
        }

        // delete silence from audio clip
        public static AudioClip TrimSilence(AudioClip clip, float min)
        {
            var samples = new float[clip.samples];

            clip.GetData(samples, 0);

            return TrimSilence(new List<float>(samples), min, clip.channels, clip.frequency);
        }

        public static AudioClip TrimSilence(List<float> samples, float min, int channels, int hz)
        {
            return TrimSilence(samples, min, channels, hz, false, false);
        }

        public static AudioClip TrimSilence(List<float> samples, float min, int channels, int hz, bool _3D, bool stream)
        {
            int i;
            //change 
            for (i = 0; i < samples.Count; i++)
            {
                if (Mathf.Abs(samples[i]) > min)
                {
                    break;
                }
            }

            samples.RemoveRange(0, i);

            for (i = samples.Count - 1; i > 0; i--)
            {
                if (Mathf.Abs(samples[i]) > min)
                {
                    break;
                }
            }

            samples.RemoveRange(i, samples.Count - i);

            var clip = AudioClip.Create("TempClip", samples.Count, channels, hz, _3D, stream);

            clip.SetData(samples.ToArray(), 0);

            return clip;
        }

        static FileStream CreateEmpty(string filepath)
        {
            var fileStream = new FileStream(filepath, FileMode.Create);
            byte emptyByte = new byte();

            for (int i = 0; i < HEADER_SIZE; i++) //preparing the header
            {
                fileStream.WriteByte(emptyByte);
            }

            return fileStream;
        }

        static void ConvertAndWrite(FileStream fileStream, AudioClip clip)
        {

            var samples = new float[clip.samples];

            clip.GetData(samples, 0);

            Int16[] intData = new Int16[samples.Length];
            //converting in 2 float[] steps to Int16[], //then Int16[] to Byte[]

            Byte[] bytesData = new Byte[samples.Length * 2];
            //bytesData array is twice the size of
            //dataSource array because a float converted in Int16 is 2 bytes.

            int rescaleFactor = 32767; //to convert float to Int16

            for (int i = 0; i < samples.Length; i++)
            {
                intData[i] = (short)(samples[i] * rescaleFactor);
                Byte[] byteArr = new Byte[2];
                byteArr = BitConverter.GetBytes(intData[i]);
                byteArr.CopyTo(bytesData, i * 2);
            }

            fileStream.Write(bytesData, 0, bytesData.Length);
        }

        static void WriteHeader(FileStream fileStream, AudioClip clip)
        {

            var hz = clip.frequency;
            var channels = clip.channels;
            var samples = clip.samples;

            fileStream.Seek(0, SeekOrigin.Begin);

            Byte[] riff = System.Text.Encoding.UTF8.GetBytes("RIFF");
            fileStream.Write(riff, 0, 4);

            Byte[] chunkSize = BitConverter.GetBytes(fileStream.Length - 8);
            fileStream.Write(chunkSize, 0, 4);

            Byte[] wave = System.Text.Encoding.UTF8.GetBytes("WAVE");
            fileStream.Write(wave, 0, 4);

            Byte[] fmt = System.Text.Encoding.UTF8.GetBytes("fmt ");
            fileStream.Write(fmt, 0, 4);

            Byte[] subChunk1 = BitConverter.GetBytes(16);
            fileStream.Write(subChunk1, 0, 4);

            UInt16 two = 2;
            UInt16 one = 1;

            Byte[] audioFormat = BitConverter.GetBytes(one);
            fileStream.Write(audioFormat, 0, 2);

            Byte[] numChannels = BitConverter.GetBytes(channels);
            fileStream.Write(numChannels, 0, 2);

            Byte[] sampleRate = BitConverter.GetBytes(hz);
            fileStream.Write(sampleRate, 0, 4);

            Byte[] byteRate = BitConverter.GetBytes(hz * channels * 2); // sampleRate * bytesPerSample*number of channels, here 44100*2*2
            fileStream.Write(byteRate, 0, 4);

            UInt16 blockAlign = (ushort)(channels * 2);
            fileStream.Write(BitConverter.GetBytes(blockAlign), 0, 2);

            UInt16 bps = 16;
            Byte[] bitsPerSample = BitConverter.GetBytes(bps);
            fileStream.Write(bitsPerSample, 0, 2);

            Byte[] datastring = System.Text.Encoding.UTF8.GetBytes("data");
            fileStream.Write(datastring, 0, 4);

            Byte[] subChunk2 = BitConverter.GetBytes(samples * channels * 2);
            fileStream.Write(subChunk2, 0, 4);

            fileStream.Close();
        }
    }

    public class Note
    {
        public String Name;
        public double Frequency;
        public int MIDINum;
        
        public Note(){}

        public Note(String name, double frequency, int midiNum)
        {
            this.Name = name;
            this.Frequency = frequency;
            this.MIDINum = midiNum;
        }

    }

    public class NoteDetector
    {
        private static readonly String[] Notes = { "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B" };
        private const float FREQUENCY_TOLERANCE = 15.0f; // cents   
        private const float A4_FREQUENCY = 440.0f;
        private const int A4_MIDI_NUMBER = 69;

        public Note GetNoteFromFrequency(float frequency){
            // confert into midi number 
            double midiNumber = A4_MIDI_NUMBER + 12 * Mathf.Log(frequency / A4_FREQUENCY)/Mathf.Log(2);
            int roundedMIDI = (int) Math.Round(midiNumber);

            // check for tolerance  
            float cents = Mathf.Abs(1200 * Mathf.Log(frequency / getFrequencyForMidiNote(roundedMIDI)));

            if(cents > FREQUENCY_TOLERANCE){
                return new Note();
            }

            return new Note(
                getNoteNameFromMidi(roundedMIDI),
                getFrequencyForMidiNote(roundedMIDI),
                roundedMIDI
            );
        }

        private float getFrequencyForMidiNote(int midiNumber){
            return A4_FREQUENCY * Mathf.Pow(2, (midiNumber - A4_MIDI_NUMBER / 12.0f));
        }
        private string getNoteNameFromMidi(int midiNumber){
            return Notes[midiNumber % 12];
        }

    }

    public class ChordDetector : MonoBehaviour
    {
        
        private readonly NoteDetector noteDetector;

        public ChordDetector()
        {
            noteDetector = new NoteDetector();
        }

        public class ChordResult
        {
            public string ChordName { get; private set; }
            public List<Note> Notes { get; private set; }
            public ChordType Type { get; private set; }

            public ChordResult(string chordName, List<Note> notes, ChordType type)
            {
                ChordName = chordName;
                Notes = notes;
                Type = type;
            }
        }

        public enum ChordType
        {
            SingleNote,
            Dyad,
            Triad,
            Seventh,
            Extended,
            Unknown
        }

        public ChordResult AnylyzeFrequencies(float[] frequencies){
            
            // check frequency list 
            if(frequencies == null || frequencies.Length == 0)
                return new ChordResult(
                    "",
                    new List<Note>(),
                    ChordType.Unknown
                );

            List<Note> detectedNotes = new List<Note>();
            foreach(float freq in frequencies)
            {
                Note note = noteDetector.GetNoteFromFrequency(freq);
                if(note.Name != null)
                {
                    detectedNotes.Add(note);
                }
            } 

            detectedNotes.OrderBy(n=> n.MIDINum).ToList();

            if(detectedNotes.Count == 0)
                return new ChordResult("", new List<Note>(), ChordType.Unknown);

            if(detectedNotes.Count == 1)
                return new ChordResult(detectedNotes[0].Name, detectedNotes, ChordType.SingleNote);

            if(detectedNotes.Count == 2)
                return AnalyzeDyad(detectedNotes);

            if(detectedNotes.Count == 3)
                return AnalyzeTriad(detectedNotes);

            
            return null;
        }

        private ChordResult AnalyzeDyad(List<Note> detectedNotes)
        {

            int interval = detectedNotes[1].MIDINum - detectedNotes[0].MIDINum;
            string intervalValue = GetIntervalValue(interval);


            return new ChordResult
            (
                intervalValue,
                detectedNotes,
                ChordType.Dyad
            );
        }
        private string GetIntervalValue(int semitones)
        {
            return semitones switch
            {
                0=> "Unison",
                1=> "Minor 2nd",
                2=> "Major 3nd",
                3=> "Minor 3rd",
                4=> "Major 3rd",
                5=> "Perfect 4th",
                6=> "Perfect 5th",
                7=> "Minor 6th",
                8=> "Major 6th",
                9=> "Minor 7th",
                10=> "Minor 7th",
                11=> "Major 7th",
                12=> "Octave",
                _=> "${semitones} semitones"

            } ;
        }

        private ChordResult AnalyzeTriad(List<Note> detectedNotes)
        {
            List<int> intervals = new List<int>();
            for(int i = 0; i < detectedNotes.Count - 1; i++ )
            {
                intervals.Add(detectedNotes[i+1].MIDINum - detectedNotes[i].MIDINum);
            }

            string chordType = DetermineTriadType(intervals);
            string rootNote = detectedNotes[0].Name;

            return new ChordResult(
                "${rootName}{chordType}",
                detectedNotes,
                ChordType.Triad
            );
        }
        private string DetermineTriadType(List<int> intervals)
    {
        if (intervals.Count >= 2)
        {
            if (intervals[0] == 4 && intervals[1] == 3) return "Major";
            if (intervals[0] == 3 && intervals[1] == 4) return "Minor";
            if (intervals[0] == 4 && intervals[1] == 4) return "Augmented";
            if (intervals[0] == 3 && intervals[1] == 3) return "Diminished";
            if (intervals[0] == 4 && intervals[1] == 2) return "sus2";
            if (intervals[0] == 3 && intervals[1] == 5) return "sus4";
        }
        return "";
    }


    }
}

