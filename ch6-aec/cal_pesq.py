from scipy.io import wavfile
import soundfile as sf
from pesq import pesq

def cal_pesq(f1,f2):
    ref,rate = sf.read(f1)
    deg,rate = sf.read(f2)
    print(f1,f2)
    print('wb', pesq(rate, ref, deg, 'wb'))
    print('nb', pesq(rate, ref, deg, 'nb'))

cal_pesq('./nearspeech.wav', './output.wav')
cal_pesq('./nearspeech.wav', './rec.wav')

