from scipy.io import wavfile
from pesq import pesq

def cal_pesq(f1,f2):
    rate, ref = wavfile.read(f1)
    rate, deg = wavfile.read(f2)
    print(f1,f2)
    print('wb', pesq(rate, ref, deg, 'wb'))
    print('nb', pesq(rate, ref, deg, 'nb'))

cal_pesq('./speech.wav', './speech_rev_1_90_4.wav')
cal_pesq('./speech.wav', './speech_rev_2_90_4.wav')
cal_pesq('./speech.wav', './speech_rev_3_90_4.wav')
cal_pesq('./speech.wav', './speech_rev_4_90_4.wav')
# cal_pesq('./speech.wav', './speech_rev_5_90_4.wav')
# cal_pesq('./speech.wav', './speech_rev_6_90_4.wav')
# cal_pesq('./speech.wav', './speech_rev_7_90_4.wav')
# cal_pesq('./speech.wav', './speech_rev_8_90_4.wav')

cal_pesq('./speech.wav', './speech_rev_ds_90_4.wav')

