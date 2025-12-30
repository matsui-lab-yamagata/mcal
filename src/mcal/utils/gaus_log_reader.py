"""Gaussian log file reader (2025/06/21)"""
import os
from pathlib import Path


def check_normal_termination(log_file: str) -> bool:
    """Check if the calculation terminated normally.

    Parameters
    ----------
    log_file : str
        Path to the log file

    Returns
    -------
    bool
        True if the calculation terminated normally, False otherwise
    """
    if not Path(log_file).exists():
        return False

    with FileReader(log_file) as f:
        line = f.reversed_readline()
        if 'Normal termination' in line:
            return True
        else:
            return False


class FileReader:
    """File reader that reads line by line from the end.

    Examples
    --------
    >>> with FileReader('sample.txt') as f:
    ...     while True:
    ...         line = f.reversed_readline()
    ...         if not line:
    ...             break
    ...         print(line)
    """
    def __init__(self, file_path, buffer_size=8192):
        self.buffer_size = buffer_size
        self.file = open(file_path, 'rb')
        # ファイルポインタをファイルの末尾に移動
        self.file.seek(0, os.SEEK_END)
        self.file_size = self.file.tell()
        # 最後のブロックの終了位置をファイルサイズに設定
        self.reversed_block_end = self.file_size
        # 最後のブロックの開始位置を計算（ファイル末尾からバッファサイズ分前、またはファイル先頭）
        self.reversed_block_start = max(0, self.reversed_block_end - self.buffer_size)
        self.lines = []

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        self.file.close()

    def reversed_readline(self):
        if not self.lines and self.reversed_block_end - self.reversed_block_start != 0:
            # ファイルポインタを現在のブロックの開始位置に設定
            self.file.seek(self.reversed_block_start)

            block = self.file.read(self.reversed_block_end - self.reversed_block_start)
            lines = block.splitlines(keepends=True)

            # 現在のブロックがファイルの最初のブロックの場合
            if self.reversed_block_start == 0:
                # ブロックの終了位置を開始位置に設定して次の読み取りを防ぐ
                self.reversed_block_end = self.reversed_block_start
            else:
                # 次のブロックの終了位置を更新(現在のブロックの最初の行の長さを加算)
                self.reversed_block_end = max(0, self.reversed_block_start + len(lines[0]))

            self.reversed_block_start = max(0, self.reversed_block_end - self.buffer_size)

            # 末尾のブロックの場合、全行を読み取る
            if self.reversed_block_end == 0:
                self.lines = lines
            else:
                self.lines = lines[1:]

            return self.lines.pop().decode('utf-8')
        elif self.lines:
            return self.lines.pop().decode('utf-8')
        else:
            return None
