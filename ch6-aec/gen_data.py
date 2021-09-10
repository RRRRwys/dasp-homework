from scipy.io import wavfile
import numpy as np

def wgn(x, snr):
    P_sig = np.sum(abs(x)**2)/len(x)
    P_noise = P_sig/10**(snr/10.0)
    return np.random.randn(len(x)) * np.sqrt(P_noise)
def norm(x):
    x = np.array(x, dtype=float)
    return (x - np.min(x))/(np.max(x)-np.min(x))

fs,speech = wavfile.read('speech.wav')
fs,speech2 = wavfile.read('speech2.wav')
fs,speech_rev = wavfile.read('speech_rev.wav')

print(speech.shape)
print(speech2.shape)
print(speech_rev.shape)

speech2 = np.append(np.zeros(speech_rev.shape[0]-speech2.shape[0]), speech2)

print(speech2.shape)

speech = norm(speech)
speech2 = norm(speech2)
speech_rev = norm(speech_rev)

rec = norm(speech2 + wgn(speech2,50) + speech_rev)
ref = speech

wavfile.write('rec.wav',fs,rec)
wavfile.write('ref.wav',fs,ref)

