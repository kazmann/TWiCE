import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# massloading.txt ファイルのパスを指定
file_path = 'massloading.txt'

# データの読み込み
data = pd.read_csv(file_path, sep='\t')  # 区切り文字がタブの場合（適宜変更）

# データの確認
data.head()

# コンター図の描画
plt.figure(figsize=(10, 8))
sns.heatmap(data.pivot_table(index='y(m)', columns='x(m)', values='ttlmassloading(kg/sq-m)'),
            cmap='viridis', annot=True, fmt='.1f', linewidths=.5)

plt.title('Mass Loading Contour')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.show()
